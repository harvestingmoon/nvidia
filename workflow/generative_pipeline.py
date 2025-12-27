"""
Generative Protein Binder Pipeline
Orchestrates the 5-step workflow using NVIDIA NIM APIs:
1. Target Structure Prediction (AlphaFold2 or OpenFold3)
2. Binder Scaffold Design (RFDiffusion)
3. Sequence Design (ProteinMPNN)
4. Complex Prediction (AlphaFold-Multimer)
5. Quality Assessment (pLDDT scoring)
"""

import os
import json
import time
import zipfile
import io
from pathlib import Path
from typing import Optional, List, Dict, Any, Tuple, Set
import requests

from .workflow_state import (
    WorkflowSession, 
    WorkflowStage, 
    StageStatus,
    TargetProteinData,
    BinderProteinData,
    ComplexAnalysisData
)


def extract_residues_from_pdb(pdb_content: str) -> Dict[str, List[int]]:
    """
    Extract residue numbers and chains from PDB content
    
    Args:
        pdb_content: PDB file content as string
        
    Returns:
        Dict mapping chain ID to list of residue numbers
    """
    residues_by_chain: Dict[str, Set[int]] = {}
    
    for line in pdb_content.split('\n'):
        if line.startswith('ATOM'):
            try:
                chain = line[21].strip() or 'A'
                res_num = int(line[22:26].strip())
                
                if chain not in residues_by_chain:
                    residues_by_chain[chain] = set()
                residues_by_chain[chain].add(res_num)
            except (ValueError, IndexError):
                continue
    
    # Convert sets to sorted lists
    return {chain: sorted(list(nums)) for chain, nums in residues_by_chain.items()}


def validate_hotspot_residues(pdb_content: str, hotspot_res: List[str]) -> Tuple[List[str], List[str]]:
    """
    Validate that hotspot residues exist in the PDB file
    
    Args:
        pdb_content: PDB file content
        hotspot_res: List of hotspot residues (e.g., ["A14", "A15"])
        
    Returns:
        Tuple of (valid_hotspots, invalid_hotspots)
    """
    residues = extract_residues_from_pdb(pdb_content)
    valid = []
    invalid = []
    
    for hotspot in hotspot_res:
        # Parse hotspot format (e.g., "A14" -> chain='A', res_num=14)
        if len(hotspot) < 2:
            invalid.append(hotspot)
            continue
            
        chain = hotspot[0].upper()
        try:
            res_num = int(hotspot[1:])
        except ValueError:
            invalid.append(hotspot)
            continue
        
        if chain in residues and res_num in residues[chain]:
            valid.append(hotspot)
        else:
            invalid.append(hotspot)
    
    return valid, invalid


def validate_and_fix_contigs(pdb_content: str, contigs: str) -> Tuple[str, List[str]]:
    """
    Validate contigs specification against PDB and auto-fix if needed.
    
    Contigs format: "A1-25/0 70-100" means:
    - A1-25: Use residues 1-25 from chain A of target
    - /0: Separator
    - 70-100: Generate 70-100 new residues for binder
    
    Args:
        pdb_content: PDB file content
        contigs: Contigs specification string
        
    Returns:
        Tuple of (fixed_contigs, list_of_warnings)
    """
    residues = extract_residues_from_pdb(pdb_content)
    warnings = []
    
    if not residues:
        return contigs, ["Could not extract residues from PDB"]
    
    # Parse contigs - format like "A1-25/0 70-100" or "A10-50/0 60-80"
    import re
    
    # Find chain-residue patterns like A1-25, B10-50
    chain_pattern = re.compile(r'([A-Z])(\d+)-(\d+)')
    
    fixed_contigs = contigs
    
    for match in chain_pattern.finditer(contigs):
        chain = match.group(1)
        start_res = int(match.group(2))
        end_res = int(match.group(3))
        original = match.group(0)
        
        if chain not in residues:
            # Chain doesn't exist, try to find any chain
            available_chain = list(residues.keys())[0] if residues else 'A'
            warnings.append(f"Chain {chain} not found in PDB, using chain {available_chain}")
            chain = available_chain
        
        if chain in residues:
            available_nums = residues[chain]
            min_res = min(available_nums)
            max_res = max(available_nums)
            
            # Check if requested range is valid
            if start_res < min_res or end_res > max_res:
                # Need to fix the range
                # Keep the same range size if possible
                range_size = end_res - start_res
                
                # Adjust to fit within available residues
                new_start = max(min_res, min(start_res, max_res - range_size))
                new_end = min(max_res, new_start + range_size)
                
                # Make sure we have at least some residues
                if new_end <= new_start:
                    new_start = min_res
                    new_end = min(max_res, min_res + range_size)
                
                new_contig = f"{chain}{new_start}-{new_end}"
                fixed_contigs = fixed_contigs.replace(original, new_contig)
                warnings.append(
                    f"Adjusted {original} to {new_contig} (PDB has {chain}{min_res}-{max_res})"
                )
    
    return fixed_contigs, warnings

def _is_zip_bytes(data: bytes) -> bool:
    return zipfile.is_zipfile(io.BytesIO(data))

def _normalize_to_pdb_text(payload) -> str:
    """
    Accepts payload as dict/list/str and returns a single PDB/CIF text string.
    Handles common shapes:
    - list[str] of PDB lines/chunks -> JOIN
    - dict with keys output_pdb/pdb/result/outputs/pdbs
    - AF3-like keys model_pdb_content/model_pdb_base64/model_cif_content
    - raw string PDB
    """
    if payload is None:
        raise ValueError("Empty payload; cannot extract PDB.")

    # list -> almost always lines/chunks; join into one text blob
    if isinstance(payload, list):
        joined = "".join(x for x in payload if isinstance(x, str))
        if not joined.strip():
            raise ValueError("List payload joined to empty text.")
        return joined

    # dict -> look for known keys
    if isinstance(payload, dict):
        if payload.get("error"):
            raise RuntimeError(f"Model error: {payload['error']}")

        # AF3-style keys (safe to support here too)
        if payload.get("model_pdb_content"):
            return payload["model_pdb_content"]
        if payload.get("model_pdb_base64"):
            return base64.b64decode(payload["model_pdb_base64"]).decode("utf-8", "ignore")
        if payload.get("model_cif_content"):
            return payload["model_cif_content"]

        # Common AF2/NIM-style keys
        for k in ("output_pdb", "pdb", "pdb_text", "result", "outputs", "pdbs"):
            if k in payload:
                v = payload[k]
                if isinstance(v, str):
                    return v
                if isinstance(v, list):
                    return "".join(x for x in v if isinstance(x, str))

        raise ValueError(f"Could not find PDB in dict payload. Keys: {list(payload.keys())}")

    # str -> could already be PDB/CIF text
    if isinstance(payload, str):
        if not payload.strip():
            raise ValueError("String payload is empty.")
        return payload

    raise TypeError(f"Unhandled payload type: {type(payload)}")

def _extract_pdb_from_zip_bytes(zbytes: bytes) -> str:
    """
    Extract PDB text from zip bytes. Supports:
    - direct *.pdb members
    - *.response member containing JSON list/dict/str
    """
    with zipfile.ZipFile(io.BytesIO(zbytes), "r") as zf:
        names = zf.namelist()

        # Prefer direct PDB file if present
        pdb_members = [n for n in names if n.lower().endswith(".pdb")]
        if pdb_members:
            raw = zf.read(pdb_members[0])
            return raw.decode("utf-8", "ignore")

        # Otherwise parse *.response JSON
        resp_members = [n for n in names if n.lower().endswith(".response")]
        if resp_members:
            raw = zf.read(resp_members[0])
            text = raw.decode("utf-8", "ignore")
            payload = json.loads(text)
            return _normalize_to_pdb_text(payload)

        raise ValueError(f"ZIP did not contain .pdb or .response members. Members: {names}")



class GenerativePipeline:
    """Orchestrates the generative protein binder design workflow"""
    
    def __init__(self, session: WorkflowSession, api_key: str, output_dir: Optional[Path] = None):
        """
        Initialize the pipeline
        
        Args:
            session: WorkflowSession to track state
            api_key: NVIDIA NGC API key
            output_dir: Directory to save output files (default: {protein_name}_{model}_output)
        """
        self.session = session
        self.api_key = api_key
        self.output_dir = output_dir or Path(f"{session.project_name}_output")
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # API endpoints
        self.endpoints = {
            "af2": "https://health.api.nvidia.com/v1/biology/deepmind/alphafold2",
            "of3": "https://health.api.nvidia.com/v1/biology/openfold/openfold3/predict",
            "rfdiffusion": "https://health.api.nvidia.com/v1/biology/ipd/rfdiffusion/generate",
            "proteinmpnn": "https://health.api.nvidia.com/v1/biology/ipd/proteinmpnn/predict",
            "multimer": "https://health.api.nvidia.com/v1/biology/deepmind/alphafold2-multimer",
            "status": "https://health.api.nvidia.com/v1/status"
        }
        
        self.headers = {
            "Authorization": f"Bearer {self.api_key}",
            "Content-Type": "application/json",
            "Accept": "application/json"
        }
    
    # ================== STEP 1: TARGET STRUCTURE PREDICTION ==================
    
    def run_target_prediction(
        self, 
        model: str = "AF2", 
        algorithm: str = "mmseqs2",
        num_diffusion_samples: int = 1,
        model_seeds: List[int] = None
    ) -> Tuple[bool, str]:
        """
        Step 1: Predict target protein 3D structure
        
        Args:
            model: "AF2" for AlphaFold2, "OF3" for OpenFold3, or "AF3" for AlphaFold3
            algorithm: MSA algorithm ("mmseqs2" or "jackhmmer") - only for AF2
            num_diffusion_samples: Number of diffusion samples for AF3 (default 1)
            model_seeds: Random seeds for AF3 reproducibility (default [42])
            
        Returns:
            (success, message) tuple
        """
        print(f"\n{'='*60}")
        print(f"STEP 1: Target Structure Prediction ({model})")
        print(f"{'='*60}")
        print(f"[DEBUG] Target sequence length: {len(self.session.target.sequence) if self.session.target.sequence else 0}")
        print(f"[DEBUG] Model: {model}, Algorithm: {algorithm}")
        
        # Validate
        if not self.session.target.sequence:
            print("[DEBUG] ERROR: No target sequence provided")
            return False, "No target sequence provided"
        
        print(f"[DEBUG] Updating stage status to IN_PROGRESS")
        self.session.update_stage_status(WorkflowStage.TARGET_PREDICTION, StageStatus.IN_PROGRESS)
        
        try:
            print(f"[DEBUG] Calling {model} API...")
            if model.upper() == "AF2":
                result = self._call_alphafold2(self.session.target.sequence, algorithm)
            elif model.upper() == "OF3":
                result = self._call_openfold3(self.session.target.sequence)
            elif model.upper() == "AF3":
                result = self._call_alphafold3(
                    self.session.target.sequence,
                    num_diffusion_samples=num_diffusion_samples,
                    model_seeds=model_seeds
                )
            else:
                print(f"[DEBUG] ERROR: Unknown model {model}")
                return False, f"Unknown model: {model}. Use 'AF2' for AlphaFold2, 'OF3' for OpenFold3, or 'AF3' for AlphaFold3."
            
            print(f"[DEBUG] API call successful, processing results...")
            
            # Store results
            if isinstance(result, list) and len(result) > 0:
                # AF2 returns list of 5 structures ranked by confidence
                print(f"[DEBUG] Received {len(result)} structures from AlphaFold2")
                self.session.target.pdb_content = result[0]  # Best (top-ranked) structure
                
                # Combine all structures with MODEL/ENDMDL markers for multi-model PDB
                all_models = []
                for idx, pdb_str in enumerate(result):
                    all_models.append(f"MODEL     {idx + 1}")
                    all_models.append(pdb_str.strip())
                    all_models.append("ENDMDL")
                self.session.target.all_structures_pdb = "\n".join(all_models)
                
                print(f"[DEBUG] Stored {len(result)} structure models")
            else:
                print(f"[DEBUG] Received single structure")
                self.session.target.pdb_content = result
            
            print(f"[DEBUG] PDB content length: {len(self.session.target.pdb_content)} characters")
            
            # Calculate pLDDT score from PDB B-factor column
            plddt_score = self._calculate_plddt(self.session.target.pdb_content)
            if plddt_score > 0:
                self.session.target.confidence_avg = plddt_score
                print(f"[DEBUG] Calculated average pLDDT: {plddt_score:.2f}")
            
            # Set model name
            model_names = {"AF2": "AlphaFold2", "OF3": "OpenFold3", "AF3": "AlphaFold3"}
            self.session.target.model_used = model_names.get(model.upper(), model)
            self.session.target.structure_predicted = True
            
            # Save to file
            filename = f"{self.session.project_name}_target_{model}.pdb"
            filepath = self.output_dir / filename
            print(f"[DEBUG] Saving to {filepath}")
            with open(filepath, 'w') as f:
                f.write(self.session.target.pdb_content)
            self.session.target.structure_file_path = str(filepath)
            
            # Save all structures if AF2 (multi-model output)
            if model.upper() == "AF2" and self.session.target.all_structures_pdb:
                all_filename = f"{self.session.project_name}_target_{model}_all.pdb"
                all_filepath = self.output_dir / all_filename
                print(f"[DEBUG] Saving all {len(result) if isinstance(result, list) else 1} structures to {all_filename}")
                with open(all_filepath, 'w') as f:
                    f.write(self.session.target.all_structures_pdb)
            
            print(f"[DEBUG] Updating stage status to COMPLETED")
            self.session.update_stage_status(WorkflowStage.TARGET_PREDICTION, StageStatus.COMPLETED)
            self.session.advance_to_stage(WorkflowStage.BINDER_SCAFFOLD_DESIGN)
            
            confidence_info = f" (pLDDT: {plddt_score:.1f})" if plddt_score > 0 else ""
            print(f"✅ Target structure predicted with {self.session.target.model_used}{confidence_info}")
            print(f"   Saved to: {filepath}")
            return True, f"Target structure predicted with {self.session.target.model_used}{confidence_info}"
            
        except Exception as e:
            print(f"[DEBUG] ERROR: {str(e)}")
            import traceback
            print(f"[DEBUG] Traceback: {traceback.format_exc()}")
            self.session.update_stage_status(WorkflowStage.TARGET_PREDICTION, StageStatus.FAILED)
            return False, f"Target prediction failed: {str(e)}"
    
    def _call_alphafold2(self, sequence: str, algorithm: str = "mmseqs2") -> List[str]:
        """
        Call AlphaFold2 API with enhanced options
        
        Args:
            sequence: Protein sequence to predict
            algorithm: MSA algorithm - "mmseqs2" (faster) or "jackhmmer" (more sensitive)
            
        Returns:
            List of 5 PDB structure strings ranked by confidence
        """
        print(f"[DEBUG] _call_alphafold2: sequence length={len(sequence)}, algorithm={algorithm}")
        print(f"[DEBUG] AlphaFold2 typically takes 5-10 minutes for structure prediction")
        
        # AlphaFold2 payload with all available options
        payload = {
            "sequence": sequence,
            "algorithm": algorithm,
            "e_value": 0.0001,          # E-value threshold for MSA
            "iterations": 1,             # Number of MSA iterations  
            "databases": ["small_bfd"],  # MSA databases to search
            "relax_prediction": False    # Whether to relax structure (slower)
        }
        
        print(f"[DEBUG] Sending POST request to {self.endpoints['af2']}")
        print(f"[DEBUG] Payload: algorithm={algorithm}, databases=small_bfd")
        
        response = requests.post(
            self.endpoints["af2"],
            headers=self.headers,
            json=payload,
            timeout=(10, 600)  # Extended timeout for long sequences
        )
        
        print(f"[DEBUG] Response status code: {response.status_code}")
        
        if response.status_code == 200:
            result = response.json()
            print(f"[DEBUG] Got immediate response, type: {type(result)}")
            # AlphaFold2 returns a list of 5 PDB strings ranked by confidence
            if isinstance(result, list):
                print(f"[DEBUG] Received {len(result)} structure predictions")
                return result
            else:
                # Handle case where result is wrapped
                return [result] if isinstance(result, str) else result
        elif response.status_code == 202:
            print(f"[DEBUG] Got 202 (Accepted), will poll for results")
            print(f"[DEBUG] AlphaFold2 is processing - this may take 5-10 minutes...")
            return self._poll_async_result(response, self.endpoints["status"])
        else:
            error_detail = response.text[:500] if len(response.text) > 500 else response.text
            print(f"[DEBUG] Error response: {error_detail}")
            raise Exception(f"AlphaFold2 API error: {response.status_code} - {error_detail}")
    
    def _call_openfold3(self, sequence: str) -> str:
        """Call OpenFold3 API"""
        msa_csv = f"key,sequence\n-1,{sequence}"
        
        payload = {
            "request_id": self.session.project_name,
            "inputs": [{
                "input_id": self.session.project_name,
                "molecules": [{
                    "type": "protein",
                    "id": "A",
                    "sequence": sequence,
                    "msa": {
                        "main_db": {
                            "csv": {
                                "alignment": msa_csv,
                                "format": "csv"
                            }
                        }
                    }
                }],
                "output_format": "pdb"
            }]
        }
        
        response = requests.post(self.endpoints["of3"], headers=self.headers, json=payload)
        
        if response.status_code == 200:
            result = response.json()
            pdb_text = (
                result.get("outputs", [{}])[0]
                .get("structures_with_scores", [{}])[0]
                .get("structure", "")
            )
            if not pdb_text.strip().startswith(("ATOM", "HETATM", "MODEL", "HEADER")):
                raise ValueError("Invalid PDB content from OpenFold3")
            return pdb_text
        else:
            raise Exception(f"OF3 API error: {response.status_code} - {response.text}")
    
    def _call_alphafold3(
        self, 
        sequence: str, 
        num_diffusion_samples: int = 1,
        model_seeds: List[int] = None
    ) -> str:
        """
        Call self-hosted AlphaFold3 server
        
        Args:
            sequence: Protein sequence to predict
            num_diffusion_samples: Number of diffusion samples (default 1)
            model_seeds: Random seeds for reproducibility (default [42])
            
        Returns:
            PDB structure string
        """
        from core.protein_models import ALPHAFOLD3_SERVER
        
        if model_seeds is None:
            model_seeds = [42]
        
        print(f"[DEBUG] _call_alphafold3: sequence length={len(sequence)}")
        print(f"[DEBUG] Server: {ALPHAFOLD3_SERVER['host']}")
        print(f"[DEBUG] Diffusion samples: {num_diffusion_samples}, Seeds: {model_seeds}")
        
        # Check server health first
        health_url = f"{ALPHAFOLD3_SERVER['host']}{ALPHAFOLD3_SERVER['health_endpoint']}"
        print(f"[DEBUG] Checking server health at {health_url}")
        
        try:
            health_response = requests.get(health_url, timeout=10)
            if health_response.status_code != 200:
                raise Exception(f"AlphaFold3 server health check failed: {health_response.status_code}")
            print(f"[DEBUG] Server health check passed")
        except requests.exceptions.RequestException as e:
            raise Exception(f"Cannot connect to AlphaFold3 server at {ALPHAFOLD3_SERVER['host']}: {str(e)}")
        
        # Prepare prediction payload
        payload = {
            "name": self.session.project_name,
            "sequences": [{"id": "A", "sequence": sequence}],
            "model_seeds": model_seeds,
            "num_diffusion_samples": num_diffusion_samples
        }
        
        predict_url = f"{ALPHAFOLD3_SERVER['host']}{ALPHAFOLD3_SERVER['predict_endpoint']}"
        print(f"[DEBUG] Sending prediction request to {predict_url}")
        
        # AlphaFold3 can take a while, use extended timeout
        response = requests.post(
            predict_url,
            headers={"Content-Type": "application/json"},
            json=payload,
            timeout=(10, 1800)  # 30 minute timeout for prediction
        )
        
        print(f"[DEBUG] Response status code: {response.status_code}")
        
        if response.status_code == 200:
            result = response.json()
            print(f"[DEBUG] Got response, keys: {result.keys() if isinstance(result, dict) else type(result)}")
            
            # Extract PDB content from response
            pdb_text = None
            
            if isinstance(result, dict):
                # Check for error first
                if result.get("error"):
                    raise Exception(f"AlphaFold3 error: {result['error']}")
                
                # Try direct PDB content keys (from our server)
                if "model_pdb_content" in result and result["model_pdb_content"]:
                    pdb_text = result["model_pdb_content"]
                    print(f"[DEBUG] Found PDB in 'model_pdb_content'")
                
                # Try base64 encoded PDB
                elif "model_pdb_base64" in result and result["model_pdb_base64"]:
                    import base64
                    pdb_text = base64.b64decode(result["model_pdb_base64"]).decode('utf-8')
                    print(f"[DEBUG] Decoded PDB from 'model_pdb_base64'")
                
                # Try CIF content and convert
                elif "model_cif_content" in result and result["model_cif_content"]:
                    cif_text = result["model_cif_content"]
                    print(f"[DEBUG] Found CIF in 'model_cif_content', using as-is (CIF format)")
                    # For now, return CIF - downstream code should handle both formats
                    pdb_text = cif_text
                
                # Try other common keys
                if not pdb_text:
                    for key in ["pdb", "structure", "output", "model", "prediction", "result"]:
                        if key in result:
                            value = result[key]
                            if isinstance(value, str) and ("ATOM" in value or "MODEL" in value):
                                pdb_text = value
                                print(f"[DEBUG] Found structure in '{key}'")
                                break
            
            elif isinstance(result, str):
                pdb_text = result
            
            if not pdb_text:
                print(f"[DEBUG] Could not parse PDB from response: {str(result)[:500]}")
                raise ValueError(f"Could not extract PDB from AlphaFold3 response. Keys: {result.keys() if isinstance(result, dict) else 'not a dict'}")
            
            # Validate content
            valid_starts = ("ATOM", "HETATM", "MODEL", "HEADER", "REMARK", "data_", "#")  # Added CIF formats
            if not pdb_text.strip().startswith(valid_starts):
                print(f"[DEBUG] Warning: Structure content may be invalid, first 100 chars: {pdb_text[:100]}")
            
            print(f"[DEBUG] Successfully extracted structure, length: {len(pdb_text)}")
            return pdb_text
            
        else:
            error_detail = response.text[:500] if len(response.text) > 500 else response.text
            print(f"[DEBUG] Error response: {error_detail}")
            raise Exception(f"AlphaFold3 API error: {response.status_code} - {error_detail}")
    
    # ================== STEP 2: BINDER SCAFFOLD DESIGN ==================
    
    def run_scaffold_design(
        self,
        contigs: str = "A1-25/0 70-100",
        hotspot_res: Optional[List[str]] = None,
        diffusion_steps: int = 15
    ) -> Tuple[bool, str]:
        """
        Step 2: Generate binder scaffold using RFDiffusion
        
        Args:
            contigs: Contig specification (e.g., "A1-25/0 70-100")
            hotspot_res: Binding hotspot residues (e.g., ["A14", "A15"])
            diffusion_steps: Number of diffusion steps
            
        Returns:
            (success, message)
        """
        print(f"\n{'='*60}")
        print(f"STEP 2: Binder Scaffold Design (RFDiffusion)")
        print(f"{'='*60}")
        print(f"[DEBUG] Contigs: {contigs}")
        print(f"[DEBUG] Hotspot residues: {hotspot_res}")
        print(f"[DEBUG] Diffusion steps: {diffusion_steps}")
        
        if not self.session.target.pdb_content:
            print(f"[DEBUG] ERROR: No target structure available")
            return False, "No target structure available. Run target prediction first."
        
        print(f"[DEBUG] Target PDB length: {len(self.session.target.pdb_content)} characters")
        print(f"[DEBUG] Updating stage status to IN_PROGRESS")
        self.session.update_stage_status(WorkflowStage.BINDER_SCAFFOLD_DESIGN, StageStatus.IN_PROGRESS)
        
        try:
            # Prepare input - use only ATOM lines, limit to first 400
            print(f"[DEBUG] Filtering ATOM lines from target PDB...")
            target_pdb_lines = [
                line for line in self.session.target.pdb_content.split("\n")
                if line.startswith("ATOM")
            ]
            print(f"[DEBUG] Found {len(target_pdb_lines)} ATOM lines, using first 400")
            target_pdb_input = "\n".join(target_pdb_lines[:400])
            
            # Get available residues for validation
            available_residues = extract_residues_from_pdb(target_pdb_input)
            print(f"[DEBUG] Available residues in PDB: {available_residues}")
            
            # Validate and fix contigs specification
            fixed_contigs, contig_warnings = validate_and_fix_contigs(target_pdb_input, contigs)
            if contig_warnings:
                for warning in contig_warnings:
                    print(f"[DEBUG] CONTIGS WARNING: {warning}")
            if fixed_contigs != contigs:
                print(f"[DEBUG] Contigs adjusted from '{contigs}' to '{fixed_contigs}'")
                contigs = fixed_contigs
            
            # Validate hotspot residues if provided
            validated_hotspots = None
            if hotspot_res:
                valid, invalid = validate_hotspot_residues(target_pdb_input, hotspot_res)
                if invalid:
                    print(f"[DEBUG] WARNING: Invalid hotspot residues removed: {invalid}")
                    # Get available residues for error message
                    available = extract_residues_from_pdb(target_pdb_input)
                    available_str = ", ".join([f"{chain}: {min(nums)}-{max(nums)}" for chain, nums in available.items()])
                    if not valid:
                        return False, f"All hotspot residues are invalid: {invalid}. Available residues: {available_str}"
                    print(f"[DEBUG] Using valid hotspots: {valid}")
                validated_hotspots = valid if valid else None
            
            # Store parameters
            self.session.binder.rfdiffusion_params = {
                "contigs": contigs,
                "hotspot_res": validated_hotspots or [],
                "diffusion_steps": diffusion_steps
            }
            
            # Call RFDiffusion
            payload = {
                "input_pdb": target_pdb_input,
                "contigs": contigs,
                "diffusion_steps": diffusion_steps
            }
            if validated_hotspots:
                payload["hotspot_res"] = validated_hotspots
            
            print(f"[DEBUG] Sending RFDiffusion request to {self.endpoints['rfdiffusion']}")
            print(f"[DEBUG] Payload: contigs={contigs}, hotspots={validated_hotspots}, steps={diffusion_steps}")
            
            response = requests.post(
                self.endpoints["rfdiffusion"],
                headers=self.headers,
                json=payload
            )
            
            print(f"[DEBUG] RFDiffusion response status: {response.status_code}")
            
            if response.status_code == 200:
                result = json.loads(response.text)
                scaffold_pdb = result["output_pdb"]
                print(f"[DEBUG] Got scaffold PDB, length: {len(scaffold_pdb)} characters")
            else:
                print(f"[DEBUG] Error response: {response.text[:200]}")
                raise Exception(f"RFDiffusion API error: {response.status_code} - {response.text}")
            
            # Store results
            self.session.binder.scaffold_pdb = scaffold_pdb
            self.session.binder.design_method = "rfdiffusion"
            
            # Save to file
            filename = f"{self.session.project_name}_scaffold_RFD.pdb"
            filepath = self.output_dir / filename
            print(f"[DEBUG] Saving scaffold to {filepath}")
            with open(filepath, 'w') as f:
                f.write(scaffold_pdb)
            self.session.binder.scaffold_file_path = str(filepath)
            
            print(f"[DEBUG] Updating stage status to COMPLETED")
            self.session.update_stage_status(WorkflowStage.BINDER_SCAFFOLD_DESIGN, StageStatus.COMPLETED)
            self.session.advance_to_stage(WorkflowStage.BINDER_SEQUENCE_DESIGN)
            
            print(f"✅ Binder scaffold generated and saved to {filepath}")
            return True, "Binder scaffold generated with RFDiffusion"
            
        except Exception as e:
            print(f"[DEBUG] ERROR: {str(e)}")
            self.session.update_stage_status(WorkflowStage.BINDER_SCAFFOLD_DESIGN, StageStatus.FAILED)
            return False, f"Scaffold design failed: {str(e)}"
    
    # ================== STEP 3: SEQUENCE DESIGN ==================
    
    def run_sequence_design(
        self,
        num_sequences: int = 10,
        sampling_temp: float = 0.1,
        use_soluble_model: bool = False
    ) -> Tuple[bool, str]:
        """
        Step 3: Design amino acid sequences for the scaffold using ProteinMPNN
        
        Args:
            num_sequences: Number of sequences to generate
            sampling_temp: Sampling temperature (lower = more conservative)
            use_soluble_model: Use soluble protein model
            
        Returns:
            (success, message)
        """
        print(f"\n{'='*60}")
        print(f"STEP 3: Sequence Design (ProteinMPNN)")
        print(f"{'='*60}")
        print(f"[DEBUG] Number of sequences: {num_sequences}")
        print(f"[DEBUG] Sampling temperature: {sampling_temp}")
        print(f"[DEBUG] Use soluble model: {use_soluble_model}")
        
        if not self.session.binder.scaffold_pdb:
            print(f"[DEBUG] ERROR: No binder scaffold available")
            return False, "No binder scaffold available. Run RFDiffusion first."
        
        print(f"[DEBUG] Scaffold PDB length: {len(self.session.binder.scaffold_pdb)} characters")
        print(f"[DEBUG] Updating stage status to IN_PROGRESS")
        self.session.update_stage_status(WorkflowStage.BINDER_SEQUENCE_DESIGN, StageStatus.IN_PROGRESS)
        
        try:
            # Prepare scaffold input
            print(f"[DEBUG] Filtering ATOM lines from scaffold PDB...")
            scaffold_lines = [
                line for line in self.session.binder.scaffold_pdb.split("\n")
                if line.startswith("ATOM")
            ]
            print(f"[DEBUG] Found {len(scaffold_lines)} ATOM lines, using first 400")
            scaffold_input = "\n".join(scaffold_lines[:400])
            
            # Call ProteinMPNN
            payload = {
                "input_pdb": scaffold_input,
                "input_pdb_chains": ["A"],
                "ca_only": False,
                "use_soluble_model": use_soluble_model,
                "num_seq_per_target": num_sequences,
                "sampling_temp": [sampling_temp]
            }
            
            print(f"[DEBUG] Sending ProteinMPNN request to {self.endpoints['proteinmpnn']}")
            print(f"[DEBUG] Payload: num_seq={num_sequences}, temp={sampling_temp}")
            
            response = requests.post(
                self.endpoints["proteinmpnn"],
                headers=self.headers,
                json=payload
            )
            
            print(f"[DEBUG] ProteinMPNN response status: {response.status_code}")
            
            if response.status_code == 200:
                result = json.loads(response.text)
                fasta_content = result["mfasta"]
                print(f"[DEBUG] Got FASTA content, length: {len(fasta_content)} characters")
            else:
                print(f"[DEBUG] Error response: {response.text[:200]}")
                raise Exception(f"ProteinMPNN API error: {response.status_code} - {response.text}")
            
            # Parse FASTA sequences
            sequences = []
            scores = []
            for line in fasta_content.split("\n"):
                if line.startswith(">"):
                    # Extract score if present
                    if "score=" in line:
                        score_str = line.split("score=")[1].split(",")[0]
                        scores.append(float(score_str))
                elif line.strip() and not line.startswith(">"):
                    sequences.append(line.strip())
            
            # Skip the first sequence (it's often the input sequence)
            if len(sequences) > 1:
                sequences = sequences[1:]
                if len(scores) > 1:
                    scores = scores[1:]
            
            # Store results
            self.session.binder.mpnn_fasta_content = fasta_content
            self.session.binder.mpnn_sequences = sequences
            self.session.binder.mpnn_scores = scores
            self.session.binder.mpnn_num_sequences = len(sequences)
            
            # Set the first sequence as default
            if sequences:
                self.session.binder.sequence = sequences[0]
                self.session.binder.selected_sequence_idx = 0
            
            # Save to file
            filename = f"{self.session.project_name}_sequences_MPNN.fa"
            filepath = self.output_dir / filename
            with open(filepath, 'w') as f:
                f.write(fasta_content)
            
            self.session.update_stage_status(WorkflowStage.BINDER_SEQUENCE_DESIGN, StageStatus.COMPLETED)
            self.session.advance_to_stage(WorkflowStage.COMPLEX_PREDICTION)
            
            print(f"✅ Generated {len(sequences)} sequences and saved to {filepath}")
            return True, f"Generated {len(sequences)} binder sequences"
            
        except Exception as e:
            self.session.update_stage_status(WorkflowStage.BINDER_SEQUENCE_DESIGN, StageStatus.FAILED)
            return False, f"Sequence design failed: {str(e)}"
    
    # ================== ALPHAFOLD3 MULTIMER ==================
    
    def _call_alphafold3_multimer(
        self,
        binder_sequence: str,
        target_sequence: str,
        num_diffusion_samples: int = 1,
        model_seeds: List[int] = None
    ) -> str:
        """
        Call AlphaFold3 server for complex (multimer) prediction
        
        Args:
            binder_sequence: Binder protein sequence
            target_sequence: Target protein sequence
            num_diffusion_samples: Number of diffusion samples
            model_seeds: Random seeds for reproducibility
            
        Returns:
            PDB structure string of the complex
        """
        from core.protein_models import ALPHAFOLD3_SERVER
        
        if model_seeds is None:
            model_seeds = [42]
        
        print(f"[DEBUG] _call_alphafold3_multimer: binder={len(binder_sequence)}AA, target={len(target_sequence)}AA")
        print(f"[DEBUG] Server: {ALPHAFOLD3_SERVER['host']}")
        
        # Check server health first
        health_url = f"{ALPHAFOLD3_SERVER['host']}{ALPHAFOLD3_SERVER['health_endpoint']}"
        print(f"[DEBUG] Checking server health at {health_url}")
        
        try:
            health_response = requests.get(health_url, timeout=10)
            if health_response.status_code != 200:
                raise Exception(f"AlphaFold3 server health check failed: {health_response.status_code}")
            print(f"[DEBUG] Server health check passed")
        except requests.exceptions.RequestException as e:
            raise Exception(f"Cannot connect to AlphaFold3 server at {ALPHAFOLD3_SERVER['host']}: {str(e)}")
        
        # Prepare multimer payload with two chains
        payload = {
            "name": f"{self.session.project_name}_complex",
            "sequences": [
                {"id": "A", "sequence": target_sequence},  # Chain A: Target
                {"id": "B", "sequence": binder_sequence}   # Chain B: Binder
            ],
            "model_seeds": model_seeds,
            "num_diffusion_samples": num_diffusion_samples
        }
        
        predict_url = f"{ALPHAFOLD3_SERVER['host']}{ALPHAFOLD3_SERVER['predict_endpoint']}"
        print(f"[DEBUG] Sending multimer prediction request to {predict_url}")
        
        # AlphaFold3 can take a while, use extended timeout
        response = requests.post(
            predict_url,
            headers={"Content-Type": "application/json"},
            json=payload,
            timeout=(10, 1800)  # 30 minute timeout for prediction
        )
        
        print(f"[DEBUG] Response status code: {response.status_code}")
        
        if response.status_code == 200:
            result = response.json()
            print(f"[DEBUG] Got response, keys: {result.keys() if isinstance(result, dict) else type(result)}")
            
            # Extract PDB content (same logic as single chain)
            pdb_text = None
            
            if isinstance(result, dict):
                if result.get("error"):
                    raise Exception(f"AlphaFold3 error: {result['error']}")
                
                # Try direct PDB content keys
                if "model_pdb_content" in result and result["model_pdb_content"]:
                    pdb_text = result["model_pdb_content"]
                    print(f"[DEBUG] Found PDB in 'model_pdb_content'")
                elif "model_pdb_base64" in result and result["model_pdb_base64"]:
                    import base64
                    pdb_text = base64.b64decode(result["model_pdb_base64"]).decode('utf-8')
                    print(f"[DEBUG] Decoded PDB from 'model_pdb_base64'")
                elif "model_cif_content" in result and result["model_cif_content"]:
                    pdb_text = result["model_cif_content"]
                    print(f"[DEBUG] Found CIF in 'model_cif_content'")
            
            if not pdb_text:
                print(f"[DEBUG] Could not parse PDB from response: {str(result)[:500]}")
                raise ValueError(f"Could not extract PDB from AlphaFold3 response. Keys: {result.keys() if isinstance(result, dict) else 'not a dict'}")
            
            print(f"[DEBUG] Successfully extracted complex structure, length: {len(pdb_text)}")
            return pdb_text
            
        else:
            error_detail = response.text[:500] if len(response.text) > 500 else response.text
            print(f"[DEBUG] Error response: {error_detail}")
            raise Exception(f"AlphaFold3 API error: {response.status_code} - {error_detail}")
    
    # ================== STEP 4: COMPLEX PREDICTION ==================
    
    def run_complex_prediction(
        self,
        sequence_idx: int = 0,
        selected_models: List[int] = [1],
        relax_prediction: bool = False,
        model_type: str = "alphafold2_multimer"  # "alphafold2_multimer" or "alphafold3"
    ) -> Tuple[bool, str]:
        """
        Step 4: Predict binder-target complex using AlphaFold-Multimer or AlphaFold3
        
        Args:
            sequence_idx: Index of the designed sequence to use
            selected_models: Which AF-Multimer models to use (1-5) - only for AF2
            relax_prediction: Whether to relax the structure - only for AF2
            model_type: "alphafold2_multimer" or "alphafold3"
            
        Returns:
            (success, message)
        """
        print(f"\n{'='*60}")
        print(f"STEP 4: Complex Prediction ({model_type})")
        print(f"{'='*60}")
        
        if not self.session.binder.mpnn_sequences:
            return False, "No designed sequences available. Run ProteinMPNN first."
        
        if sequence_idx >= len(self.session.binder.mpnn_sequences):
            return False, f"Invalid sequence index: {sequence_idx}"
        
        self.session.update_stage_status(WorkflowStage.COMPLEX_PREDICTION, StageStatus.IN_PROGRESS)
        
        try:
            # Get the selected binder sequence
            binder_seq = self.session.binder.mpnn_sequences[sequence_idx]
            target_seq = self.session.target.sequence
            
            self.session.binder.selected_sequence_idx = sequence_idx
            self.session.binder.sequence = binder_seq
            
            complex_pdb = None
            
            if model_type == "alphafold3":
                # Use AlphaFold3 server for complex prediction
                print(f"[DEBUG] Using AlphaFold3 for complex prediction")
                complex_pdb = self._call_alphafold3_multimer(
                    binder_sequence=binder_seq,
                    target_sequence=target_seq,
                    num_diffusion_samples=1,
                    model_seeds=[42]
                )
                self.session.complex.docking_method = "alphafold3"
                
            else:
                # Use AlphaFold2-Multimer (NVIDIA API)
                print(f"[DEBUG] Using AlphaFold2-Multimer for complex prediction")
                payload = {
                    "sequences": [binder_seq, target_seq],
                    "selected_models": selected_models,
                    "relax_prediction": relax_prediction,
                    "databases": ["small_bfd"]
                }
                
                response = requests.post(
                    self.endpoints["multimer"],
                    headers=self.headers,
                    json=payload,
                    timeout=(10, 310)
                )
                
                if response.status_code == 200:
                    result = response
                elif response.status_code == 202:
                    result = self._poll_async_multimer(response)
                else:
                    raise Exception(f"Multimer API error: {response.status_code} - {response.text}")

                content = result.content

                raw_filename = f"{self.session.project_name}_multimer_{sequence_idx+1}.bin"
                raw_path = self.output_dir / raw_filename
                with open(raw_path, "wb") as f:
                    f.write(content)

                # If Zip file
                if _is_zip_bytes(content):
                    zip_filename = f"{self.session.project_name}_multimer_{sequence_idx+1}.zip"
                    zip_path = self.output_dir / zip_filename
                    with open(zip_path, "wb") as f:
                        f.write(content)
                    self.session.complex.multimer_zip_path = str(zip_path)

                    complex_pdb = _extract_pdb_from_zip_bytes(content)

                else:
                    # Not a zip: JSON or raw text
                    try:
                        payload = result.json()
                    except Exception:
                        # fall back to text (decode bytes safely)
                        payload = result.text

                    complex_pdb = _normalize_to_pdb_text(payload)
                    
                # # Save ZIP file
                # zip_filename = f"{self.session.project_name}_multimer_{sequence_idx+1}.zip"
                # zip_path = self.output_dir / zip_filename
                # with open(zip_path, 'wb') as f:
                #     f.write(result.content)
                # self.session.complex.multimer_zip_path = str(zip_path)
                
                # # Extract PDB from ZIP
                # with zipfile.ZipFile(zip_path, 'r') as zf:
                #     names = zf.namelist()
                #     for name in names:
                #         raw = zf.read(name)
                #         text = raw.decode("utf-8", "ignore")
                #         obj = json.loads(text)
                #         complex_pdb = "".join(obj) if isinstance(obj, list) else str(obj)
                
                # self.session.complex.docking_method = "alphafold2_multimer"
                # self.session.complex.multimer_model_used = selected_models[0] if selected_models else 1
            
            # Store results
            self.session.complex.complex_pdb = complex_pdb
            self.session.complex.docking_method = "alphafold2_multimer"
            self.session.complex.multimer_model_used = selected_models[0] if selected_models else 1
            
            # Save PDB
            pdb_filename = f"{self.session.project_name}_complex_{sequence_idx+1}.pdb"
            pdb_path = self.output_dir / pdb_filename
            with open(pdb_path, 'w') as f:
                f.write(complex_pdb)
            
            # Calculate pLDDT score
            plddt = self._calculate_plddt(complex_pdb)
            self.session.complex.plddt_score = plddt
            self.session.complex.quality_score = int(plddt)
            
            # Assign grade
            if plddt > 90:
                self.session.complex.quality_grade = "Excellent"
            elif plddt > 70:
                self.session.complex.quality_grade = "Good"
            elif plddt > 50:
                self.session.complex.quality_grade = "Fair"
            else:
                self.session.complex.quality_grade = "Poor"
            
            self.session.update_stage_status(WorkflowStage.COMPLEX_PREDICTION, StageStatus.COMPLETED)
            self.session.advance_to_stage(WorkflowStage.RESULTS)
            
            print(f"✅ Complex predicted with pLDDT={plddt:.2f} and saved to {pdb_path}")
            return True, f"Complex predicted with pLDDT score: {plddt:.2f}"
            
        except Exception as e:
            self.session.update_stage_status(WorkflowStage.COMPLEX_PREDICTION, StageStatus.FAILED)
            return False, f"Complex prediction failed: {str(e)}"
    
    # ================== STEP 5: BATCH COMPLEX PREDICTION ==================
    
    def run_batch_complex_prediction(
        self,
        num_candidates: int = 3,
        selected_models: List[int] = [1],
        relax_prediction: bool = False,
        model_type: str = "alphafold2_multimer"  # "alphafold2_multimer" or "alphafold3"
    ) -> Tuple[bool, str]:
        """
        Step 5: Predict complexes for multiple sequence candidates and rank them
        
        Args:
            num_candidates: Number of top sequences to evaluate
            selected_models: Which AF-Multimer models to use
            relax_prediction: Whether to relax structures
            model_type: "alphafold2_multimer" or "alphafold3"
            
        Returns:
            (success, message)
        """
        print(f"\n{'='*60}")
        print(f"STEP 5: Batch Complex Prediction & Ranking ({model_type})")
        print(f"{'='*60}")
        
        num_candidates = min(num_candidates, len(self.session.binder.mpnn_sequences))
        rankings = []
        
        for i in range(num_candidates):
            print(f"\nEvaluating candidate {i+1}/{num_candidates}...")
            success, msg = self.run_complex_prediction(
                i, selected_models, relax_prediction, model_type=model_type
            )
            
            if success:
                rankings.append({
                    "sequence_idx": i,
                    "sequence": self.session.binder.mpnn_sequences[i],
                    "plddt_score": self.session.complex.plddt_score,
                    "quality_grade": self.session.complex.quality_grade,
                    "pdb_path": str(self.output_dir / f"{self.session.project_name}_complex_{i+1}.pdb")
                })
        
        # Sort by pLDDT score
        rankings.sort(key=lambda x: x["plddt_score"], reverse=True)
        self.session.complex.candidate_rankings = rankings
        
        # Save rankings to file
        rankings_file = self.output_dir / "pLDDT_scores.txt"
        with open(rankings_file, 'w') as f:
            for rank in rankings:
                f.write(f"{Path(rank['pdb_path']).name}\t{rank['plddt_score']:.6f}\n")
        
        print(f"\n✅ Evaluated {len(rankings)} candidates. Best pLDDT: {rankings[0]['plddt_score']:.2f}")
        return True, f"Evaluated {len(rankings)} candidates"
    
    # ================== HELPER METHODS ==================
    
    def _poll_async_result(self, initial_response: requests.Response, status_endpoint: str) -> Any:
        """
        Poll for async API results with detailed progress tracking
        
        Args:
            initial_response: The initial 202 response from the API
            status_endpoint: The endpoint to poll for status
            
        Returns:
            The result JSON from the API
        """
        req_id = initial_response.headers.get("NVCF-REQID") or initial_response.headers.get("nvcf-reqid")
        if not req_id:
            raise Exception("No request ID in async response")
        
        print(f"\n{'='*50}")
        print(f"⏳ ASYNC REQUEST SUBMITTED")
        print(f"   Request ID: {req_id}")
        print(f"   Status Endpoint: {status_endpoint}")
        print(f"{'='*50}")
        
        max_attempts = 180  # 30 minutes at 10-second intervals
        poll_interval = 10
        
        # Progress stages for AlphaFold2
        stages = [
            (0, 30, "🔍 Running MSA search..."),
            (30, 60, "📊 Building multiple sequence alignment..."),
            (60, 120, "🧬 Running neural network structure prediction..."),
            (120, 180, "⚙️ Refining structure predictions...")
        ]
        
        for attempt in range(max_attempts):
            elapsed_seconds = (attempt + 1) * poll_interval
            elapsed_minutes = elapsed_seconds // 60
            elapsed_secs = elapsed_seconds % 60
            
            # Determine current stage message
            stage_msg = "Processing..."
            for start, end, msg in stages:
                if start <= elapsed_seconds // 10 < end:
                    stage_msg = msg
                    break
            
            time.sleep(poll_interval)
            
            try:
                result = requests.get(
                    f"{status_endpoint}/{req_id}",
                    headers=self.headers,
                    timeout=(10, 310)
                )
                
                if result.status_code == 200:
                    print(f"\n✅ Request completed after {elapsed_minutes}m {elapsed_secs}s")
                    print(f"   Total polling attempts: {attempt + 1}")
                    return result.json()
                    
                elif result.status_code == 202:
                    # Print progress every 30 seconds
                    if (attempt + 1) % 3 == 0:
                        progress_pct = min(95, int((elapsed_seconds / 600) * 100))  # Cap at 95%
                        print(f"   [{elapsed_minutes:02d}:{elapsed_secs:02d}] {stage_msg} ({progress_pct}%)")
                    continue
                    
                else:
                    error_text = result.text[:300] if len(result.text) > 300 else result.text
                    raise Exception(f"Polling failed: {result.status_code} - {error_text}")
                    
            except requests.exceptions.Timeout:
                print(f"   [WARN] Poll request timed out at attempt {attempt + 1}, retrying...")
                continue
            except requests.exceptions.RequestException as e:
                if attempt < max_attempts - 1:
                    print(f"   [WARN] Request error: {str(e)}, retrying...")
                    continue
                raise e
        
        timeout_minutes = max_attempts * poll_interval // 60
        raise Exception(
            f"Request timed out after {timeout_minutes} minutes.\n"
            f"This can happen with:\n"
            f"• Very long sequences (>500 residues)\n"
            f"• Server under heavy load\n"
            f"• Complex protein structures\n\n"
            f"💡 Try: shorter sequence, OpenFold3, or try again later."
        )
    
    def _poll_async_multimer(self, initial_response: requests.Response) -> requests.Response:
        """Poll for AlphaFold-Multimer results (returns raw response with ZIP)"""
        req_id = initial_response.headers.get("NVCF-REQID") or initial_response.headers.get("nvcf-reqid")
        if not req_id:
            raise Exception("No request ID in multimer response")
        
        print(f"⏳ Polling for multimer results (request ID: {req_id})...")
        
        max_attempts = 180
        poll_interval = 10
        
        for attempt in range(max_attempts):
            time.sleep(poll_interval)
            
            try:
                result = requests.get(
                    f"{self.endpoints['status']}/{req_id}",
                    headers=self.headers,
                    timeout=(10, 310)
                )
                
                if result.status_code == 200:
                    print(f"✅ Multimer completed after {(attempt + 1) * poll_interval}s")
                    return result
                elif result.status_code == 202:
                    if (attempt + 1) % 6 == 0:
                        print(f"   Still polling... ({(attempt + 1) * poll_interval}s)")
                    continue
                else:
                    raise Exception(f"Polling failed: {result.status_code} - {result.text}")
                    
            except requests.exceptions.RequestException as e:
                if attempt < max_attempts - 1:
                    continue
                raise e
        
        raise Exception(f"Timeout after {max_attempts * poll_interval // 60} minutes")
    
    def _calculate_plddt(self, pdb_string: str) -> float:
        """Calculate average pLDDT score from PDB B-factor column"""
        total_plddt = 0.0
        atom_count = 0
        
        for line in pdb_string.splitlines():
            if line.startswith("ATOM"):
                atom_name = line[12:16].strip()
                if atom_name == "CA":  # Only alpha carbons
                    try:
                        plddt = float(line[60:66].strip())
                        total_plddt += plddt
                        atom_count += 1
                    except (ValueError, IndexError):
                        continue
        
        return total_plddt / atom_count if atom_count > 0 else 0.0
    
    # ================== CONVENIENCE METHODS ==================
    
    def run_full_pipeline(
        self,
        model: str = "AF2",
        contigs: str = "A1-25/0 70-100",
        hotspot_res: Optional[List[str]] = None,
        num_sequences: int = 10,
        num_candidates: int = 1
    ) -> Tuple[bool, str]:
        """
        Run the complete pipeline from start to finish
        
        Args:
            model: "AF2" or "OF3" for target prediction
            contigs: RFDiffusion contig specification
            hotspot_res: Binding hotspot residues
            num_sequences: Number of sequences to design
            num_candidates: Number of top candidates to evaluate
            
        Returns:
            (success, message)
        """
        print(f"\n{'#'*60}")
        print(f"# GENERATIVE PROTEIN BINDER DESIGN PIPELINE")
        print(f"# Project: {self.session.project_name}")
        print(f"{'#'*60}\n")
        
        # Step 1: Target prediction
        success, msg = self.run_target_prediction(model)
        if not success:
            return False, f"Pipeline failed at Step 1: {msg}"
        
        # Step 2: Scaffold design
        success, msg = self.run_scaffold_design(contigs, hotspot_res)
        if not success:
            return False, f"Pipeline failed at Step 2: {msg}"
        
        # Step 3: Sequence design
        success, msg = self.run_sequence_design(num_sequences)
        if not success:
            return False, f"Pipeline failed at Step 3: {msg}"
        
        # Step 4 & 5: Complex prediction (batch if num_candidates > 1)
        if num_candidates > 1:
            success, msg = self.run_batch_complex_prediction(num_candidates)
        else:
            success, msg = self.run_complex_prediction(0)
        
        if not success:
            return False, f"Pipeline failed at Step 4/5: {msg}"
        
        print(f"\n{'#'*60}")
        print(f"# PIPELINE COMPLETE!")
        print(f"# Output directory: {self.output_dir}")
        print(f"{'#'*60}\n")
        
        return True, "Full pipeline completed successfully"
    
    def get_stage_summary(self) -> Dict[str, Any]:
        """Get summary of current pipeline state"""
        return {
            "current_stage": self.session.current_stage.value,
            "stage_statuses": self.session.stage_statuses,
            "target_predicted": self.session.target.structure_predicted,
            "scaffold_generated": bool(self.session.binder.scaffold_pdb),
            "sequences_designed": len(self.session.binder.mpnn_sequences),
            "complex_predicted": bool(self.session.complex.complex_pdb),
            "plddt_score": self.session.complex.plddt_score,
            "output_directory": str(self.output_dir)
        }
