"""
Protein Structure Prediction with NVIDIA AlphaFold2
A Streamlit application for predicting protein 3D structures using NVIDIA NIM (AlphaFold2 model)
"""

import streamlit as st
import requests
import json
import os
import tempfile
from typing import Optional, Dict, Any
import py3Dmol
import streamlit.components.v1 as components
from langchain.llms.base import LLM
from langchain.chains import LLMChain
from langchain.prompts import PromptTemplate
from langchain.schema import Generation, LLMResult
from pydantic import Field
import re
import time

# Custom NVIDIA NIM LangChain LLM wrapper
class NVIDIANIMLangChainLLM(LLM):
    """
    Custom LangChain LLM wrapper for NVIDIA NIM API
    """
    
    nim_endpoint: str = Field(...)
    api_key: str = Field(...)
    model_name: str = "e3dfc6dd-fc27-4f0e-9ede-94412256af18"  # ai-alphafold2 function ID
    
    @property
    def _llm_type(self) -> str:
        return "nvidia_nim"
    
    def _call(self, prompt: str, stop: Optional[list] = None) -> str:
        """
        Call the NVIDIA NIM API with the protein sequence
        """
        headers = {
            "Authorization": f"Bearer {self.api_key}",
            "Content-Type": "application/json"
        }
        
        # Prepare the payload for alphafold2 model
        payload = {
            "sequence": prompt.strip()
        }
        
        try:
            # Use the NVIDIA Cloud Functions API format
            response = requests.post(
                f"{self.nim_endpoint}/v2/nvcf/pexec/functions/{self.model_name}",
                headers=headers,
                json=payload,
                timeout=300  # 5 minute timeout for structure prediction
            )
            
            if response.status_code == 200:
                result = response.json()
                # Extract PDB content from the response
                if "pdb" in result:
                    return result["pdb"]
                elif "output" in result and "pdb" in result["output"]:
                    return result["output"]["pdb"]
                else:
                    return str(result)  # Return raw response if PDB not found in expected location
            elif response.status_code == 202:
                # Handle asynchronous processing
                result = response.json()
                if "reqId" in result:
                    return self._poll_for_result(result["reqId"])
                else:
                    return f"Task submitted successfully but no request ID received: {response.text}"
            else:
                st.error(f"API Error: {response.status_code} - {response.text}")
                return f"Error: {response.status_code} - {response.text}"
                
        except requests.exceptions.RequestException as e:
            st.error(f"Request failed: {str(e)}")
            return f"Request failed: {str(e)}"
    
    def _poll_for_result(self, request_id: str, max_attempts: int = 60) -> str:
        """
        Poll for the result of an asynchronous request
        
        Args:
            request_id (str): The request ID to poll for
            max_attempts (int): Maximum number of polling attempts
            
        Returns:
            str: The final result or error message
        """
        import time
        
        headers = {
            "Authorization": f"Bearer {self.api_key}",
            "Content-Type": "application/json"
        }
        
        for attempt in range(max_attempts):
            try:
                # Poll the status endpoint
                poll_response = requests.get(
                    f"{self.nim_endpoint}/v2/nvcf/pexec/status/{request_id}",
                    headers=headers,
                    timeout=30
                )
                
                if poll_response.status_code == 200:
                    result = poll_response.json()
                    status = result.get("status", "").upper()
                    
                    if status == "COMPLETED":
                        # Get the final result
                        result_response = requests.get(
                            f"{self.nim_endpoint}/v2/nvcf/pexec/response/{request_id}",
                            headers=headers,
                            timeout=30
                        )
                        
                        if result_response.status_code == 200:
                            final_result = result_response.json()
                            # Extract PDB content from the final result
                            if "pdb" in final_result:
                                return final_result["pdb"]
                            elif "output" in final_result:
                                if isinstance(final_result["output"], dict) and "pdb" in final_result["output"]:
                                    return final_result["output"]["pdb"]
                                elif isinstance(final_result["output"], str):
                                    return final_result["output"]
                            return str(final_result)
                        else:
                            return f"Failed to get result: {result_response.status_code} - {result_response.text}"
                    
                    elif status == "FAILED":
                        error_msg = result.get("error", "Unknown error")
                        return f"Prediction failed: {error_msg}"
                    
                    elif status in ["PENDING", "IN_PROGRESS", "QUEUED"]:
                        # Update progress in Streamlit
                        progress_msg = f"Status: {status} (attempt {attempt + 1}/{max_attempts})"
                        st.info(progress_msg)
                        time.sleep(10)  # Wait 10 seconds before next poll
                        continue
                    
                    else:
                        return f"Unknown status: {status}"
                
                else:
                    time.sleep(10)  # Wait and retry on non-200 status
                    continue
                    
            except requests.exceptions.RequestException as e:
                if attempt == max_attempts - 1:
                    return f"Polling failed after {max_attempts} attempts: {str(e)}"
                time.sleep(10)
                continue
        
        return f"Timeout: Prediction did not complete within {max_attempts * 10} seconds"
    
    def _generate(self, prompts: list, stop: Optional[list] = None) -> LLMResult:
        """Generate method for LangChain compatibility"""
        generations = []
        for prompt in prompts:
            result = self._call(prompt, stop)
            generations.append([Generation(text=result)])
        return LLMResult(generations=generations)

def validate_protein_sequence(sequence: str) -> tuple[bool, str]:
    """
    Validate if the input is a valid amino acid sequence
    
    Args:
        sequence (str): Input protein sequence
        
    Returns:
        tuple: (is_valid, error_message)
    """
    if not sequence or len(sequence.strip()) == 0:
        return False, "Please enter a protein sequence"
    
    # Remove whitespace and convert to uppercase
    clean_sequence = re.sub(r'\s+', '', sequence.upper())
    
    # Check if sequence contains only valid amino acid codes
    valid_amino_acids = set('ACDEFGHIKLMNPQRSTVWY')
    invalid_chars = set(clean_sequence) - valid_amino_acids
    
    if invalid_chars:
        return False, f"Invalid amino acid characters found: {', '.join(invalid_chars)}"
    
    # Check minimum length (typically at least 10 amino acids for meaningful structure prediction)
    if len(clean_sequence) < 10:
        return False, "Sequence too short. Please enter at least 10 amino acids"
    
    # Check maximum length (AlphaFold2 works well with sequences up to ~2000 residues)
    if len(clean_sequence) > 2000:
        return False, "Sequence too long. Please enter less than 2000 amino acids"
    
    return True, clean_sequence

def create_3d_visualization(pdb_content: str) -> str:
    """
    Create a 3D molecular visualization using py3Dmol
    
    Args:
        pdb_content (str): PDB file content as string
        
    Returns:
        str: HTML string for the 3D visualization
    """
    # Create py3Dmol viewer
    viewer = py3Dmol.view(width=800, height=600)
    viewer.addModel(pdb_content, 'pdb')
    viewer.setStyle({'cartoon': {'color': 'spectrum'}})
    viewer.zoomTo()
    viewer.spin(True)
    
    # Return the HTML representation
    return viewer._make_html()

def save_pdb_file(pdb_content: str, filename: str = "predicted_structure.pdb") -> str:
    """
    Save PDB content to a temporary file
    
    Args:
        pdb_content (str): PDB file content
        filename (str): Filename for the PDB file
        
    Returns:
        str: Path to the saved file
    """
    temp_dir = tempfile.gettempdir()
    filepath = os.path.join(temp_dir, filename)
    
    with open(filepath, 'w') as f:
        f.write(pdb_content)
    
    return filepath

def setup_langchain_pipeline(nim_endpoint: str, api_key: str) -> LLMChain:
    """
    Set up the LangChain pipeline with NVIDIA NIM
    
    Args:
        nim_endpoint (str): NVIDIA NIM endpoint URL
        api_key (str): API key for authentication
        
    Returns:
        LLMChain: Configured LangChain pipeline
    """
    # Create the custom LLM
    llm = NVIDIANIMLangChainLLM(
        nim_endpoint=nim_endpoint,
        api_key=api_key
    )
    
    # Create a simple prompt template
    prompt_template = PromptTemplate(
        input_variables=["sequence"],
        template="{sequence}"
    )
    
    # Create the chain
    chain = LLMChain(llm=llm, prompt=prompt_template)
    
    return chain

def main():
    """
    Main Streamlit application
    """
    st.set_page_config(
        page_title="Protein Structure Prediction with NVIDIA AlphaFold2",
        page_icon="🧬",
        layout="wide"
    )
    
    # Application title and description
    st.title("🧬 Protein Structure Prediction with NVIDIA AlphaFold2")
    st.markdown("""
    This application uses NVIDIA's AlphaFold2 model to predict protein 3D structures from amino acid sequences.
    Simply paste your protein sequence below and click 'Predict Structure' to generate the 3D model.
    """)
    
    # Sidebar for configuration
    st.sidebar.header("⚙️ Configuration")
    
    # Environment variables for NIM endpoint and API key
    default_endpoint = os.getenv("NIM_ENDPOINT", "https://api.nvcf.nvidia.com")
    default_api_key = os.getenv("NVIDIA_API_KEY", "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")
    
    nim_endpoint = st.sidebar.text_input(
        "NVIDIA NIM Endpoint",
        value=default_endpoint,
        help="The NVIDIA NIM API endpoint URL"
    )
    
    api_key = st.sidebar.text_input(
        "API Key",
        value=default_api_key,
        type="password",
        help="Your NVIDIA API key"
    )
    
    # Main interface
    col1, col2 = st.columns([1, 1])
    
    with col1:
        st.header("📝 Input")
        
        # Text area for protein sequence input
        sequence_input = st.text_area(
            "Enter Protein Amino Acid Sequence:",
            height=200,
            placeholder="Example: MDSKGSSQKGSRLLLLLVVSNLLLCQGVVSTPKDYFVTKELVDRIRDCLNFAKTGNTVEKRQVLNLEKEMDLKQIQDDLGFTDYQEIKNTGLYKQWLLLRRKFIQRLQLRVTNDVFHSYGRCLKKLAAYQGPYGNIVYLASILQTQKFIMLLITQKIYAKNHSLKVMKLWQDDMKLYFLGNLKTGHDTNVLTAFKDKDYNQSSNFYNMHQYPPHDLKITLLRKKKWVDIVDSIELVTLEQAKKYYQYFKQKPSKMFLYGLQFKEMVPEIMSTSELLALRQFGLTYLKWKYVLKDEGKQFYQKTVFKDPDDQTLKYLVLFADYVWQVNYCSSITKYNLKQQQQQYDLTIAQFVTLAKLKQNLAMLQITKKLTMKPNKVPEIIDNVLLMNSAGDLEPDYFVTKLVD"
        )
        
        # Example sequences
        st.subheader("📚 Example Sequences")
        example_sequences = {
            "Short Example (Insulin B-chain)": "FVNQHLCGSHLVEALYLVCGERGFFYTPKT",
            "Medium Example (Lysozyme fragment)": "KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL",
            "Sample Sequence": "MDSKGSSQKGSRLLLLLVVSNLLLCQGVVST"
        }
        
        selected_example = st.selectbox("Choose an example:", list(example_sequences.keys()))
        if st.button("Load Example"):
            st.rerun()
        
        # Use example if selected
        if selected_example and st.session_state.get('load_example'):
            sequence_input = example_sequences[selected_example]
        
        # Predict button
        predict_button = st.button("🔬 Predict Structure", type="primary")
    
    with col2:
        st.header("📊 Results")
        
        if predict_button:
            if not nim_endpoint or not api_key:
                st.error("Please provide both NIM endpoint and API key in the sidebar.")
                return
            
            # Validate the input sequence
            is_valid, result = validate_protein_sequence(sequence_input)
            
            if not is_valid:
                st.error(f"Invalid sequence: {result}")
                return
            
            clean_sequence = result
            st.success(f"Valid protein sequence with {len(clean_sequence)} amino acids")
            
            # Show progress
            with st.spinner("Predicting protein structure... This may take a few minutes."):
                try:
                    # Set up LangChain pipeline
                    chain = setup_langchain_pipeline(nim_endpoint, api_key)
                    
                    # Run prediction
                    pdb_result = chain.run(sequence=clean_sequence)
                    
                    # Check if the result looks like a PDB file
                    if pdb_result.startswith("Error:") or "Request failed" in pdb_result:
                        st.error(f"Prediction failed: {pdb_result}")
                        return
                    
                    # Store results in session state
                    st.session_state['pdb_content'] = pdb_result
                    st.session_state['sequence'] = clean_sequence
                    
                    st.success("Structure prediction completed!")
                    
                except Exception as e:
                    st.error(f"An error occurred during prediction: {str(e)}")
                    return
    
    # Display results if available
    if 'pdb_content' in st.session_state:
        st.header("🧬 3D Structure Visualization")
        
        # Create two columns for visualization and download
        viz_col1, viz_col2 = st.columns([3, 1])
        
        with viz_col1:
            try:
                # Create 3D visualization
                html_content = create_3d_visualization(st.session_state['pdb_content'])
                components.html(html_content, height=600)
            except Exception as e:
                st.error(f"Visualization error: {str(e)}")
                st.text_area("Raw PDB Content:", st.session_state['pdb_content'], height=300)
        
        with viz_col2:
            st.subheader("📥 Download")
            
            # Save PDB file
            try:
                filename = f"protein_structure_{len(st.session_state.get('sequence', ''))}_aa.pdb"
                
                # Provide download button
                st.download_button(
                    label="Download PDB File",
                    data=st.session_state['pdb_content'],
                    file_name=filename,
                    mime="chemical/x-pdb"
                )
                
                st.info(f"Sequence length: {len(st.session_state.get('sequence', ''))} amino acids")
                
            except Exception as e:
                st.error(f"Download preparation failed: {str(e)}")
        
        # Show PDB content preview
        with st.expander("📄 View Raw PDB Content"):
            st.text_area("PDB File Content:", st.session_state['pdb_content'], height=300)

if __name__ == "__main__":
    main()
