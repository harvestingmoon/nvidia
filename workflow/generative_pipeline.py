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
from typing import Optional, List, Dict, Any, Tuple
import requests

from .workflow_state import (
    WorkflowSession, 
    WorkflowStage, 
    StageStatus,
    TargetProteinData,
    BinderProteinData,
    ComplexAnalysisData
)


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
    
    def run_target_prediction(self, model: str = "AF2", algorithm: str = "mmseqs2") -> Tuple[bool, str]:
        """
        Step 1: Predict target protein 3D structure
        
        Args:
            model: "AF2" for AlphaFold2 or "OF3" for OpenFold3
            algorithm: MSA algorithm ("mmseqs2" or "jackhmmer")
            
        Returns:
            (success, message)
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
            else:
                print(f"[DEBUG] ERROR: Unknown model {model}")
                return False, f"Unknown model: {model}"
            
            print(f"[DEBUG] API call successful, processing results...")
            
            # Store results
            if isinstance(result, list) and len(result) > 0:
                # AF2 returns list of 5 structures
                print(f"[DEBUG] Received {len(result)} structures from AF2")
                self.session.target.pdb_content = result[0]
                self.session.target.all_structures_pdb = "\n".join(result)
            else:
                print(f"[DEBUG] Received single structure")
                self.session.target.pdb_content = result
            
            print(f"[DEBUG] PDB content length: {len(self.session.target.pdb_content)} characters")
            
            self.session.target.model_used = model.upper()
            self.session.target.structure_predicted = True
            
            # Save to file
            filename = f"{self.session.project_name}_target_{model}.pdb"
            filepath = self.output_dir / filename
            print(f"[DEBUG] Saving to {filepath}")
            with open(filepath, 'w') as f:
                f.write(self.session.target.pdb_content)
            self.session.target.structure_file_path = str(filepath)
            
            # Save all structures if AF2
            if model.upper() == "AF2" and self.session.target.all_structures_pdb:
                all_filename = f"{self.session.project_name}_target_{model}_all.pdb"
                print(f"[DEBUG] Saving all 5 structures to {all_filename}")
                with open(self.output_dir / all_filename, 'w') as f:
                    f.write(self.session.target.all_structures_pdb)
            
            print(f"[DEBUG] Updating stage status to COMPLETED")
            self.session.update_stage_status(WorkflowStage.TARGET_PREDICTION, StageStatus.COMPLETED)
            self.session.advance_to_stage(WorkflowStage.BINDER_SCAFFOLD_DESIGN)
            
            print(f"✅ Target structure predicted and saved to {filepath}")
            return True, f"Target structure predicted with {model}"
            
        except Exception as e:
            print(f"[DEBUG] ERROR: {str(e)}")
            self.session.update_stage_status(WorkflowStage.TARGET_PREDICTION, StageStatus.FAILED)
            return False, f"Target prediction failed: {str(e)}"
    
    def _call_alphafold2(self, sequence: str, algorithm: str = "mmseqs2") -> List[str]:
        """Call AlphaFold2 API"""
        print(f"[DEBUG] _call_alphafold2: sequence length={len(sequence)}, algorithm={algorithm}")
        
        payload = {
            "sequence": sequence,
            "algorithm": algorithm
        }
        
        print(f"[DEBUG] Sending POST request to {self.endpoints['af2']}")
        response = requests.post(
            self.endpoints["af2"],
            headers=self.headers,
            json=payload,
            timeout=(10, 310)
        )
        
        print(f"[DEBUG] Response status code: {response.status_code}")
        
        if response.status_code == 200:
            result = response.json()
            print(f"[DEBUG] Got immediate response, type: {type(result)}")
            return result  # List of 5 PDB strings
        elif response.status_code == 202:
            print(f"[DEBUG] Got 202, will poll for results")
            return self._poll_async_result(response, self.endpoints["status"])
        else:
            print(f"[DEBUG] Error response: {response.text[:200]}")
            raise Exception(f"AF2 API error: {response.status_code} - {response.text}")
    
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
            
            # Store parameters
            self.session.binder.rfdiffusion_params = {
                "contigs": contigs,
                "hotspot_res": hotspot_res or [],
                "diffusion_steps": diffusion_steps
            }
            
            # Call RFDiffusion
            payload = {
                "input_pdb": target_pdb_input,
                "contigs": contigs,
                "diffusion_steps": diffusion_steps
            }
            if hotspot_res:
                payload["hotspot_res"] = hotspot_res
            
            print(f"[DEBUG] Sending RFDiffusion request to {self.endpoints['rfdiffusion']}")
            print(f"[DEBUG] Payload: contigs={contigs}, hotspots={hotspot_res}, steps={diffusion_steps}")
            
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
    
    # ================== STEP 4: COMPLEX PREDICTION ==================
    
    def run_complex_prediction(
        self,
        sequence_idx: int = 0,
        selected_models: List[int] = [1],
        relax_prediction: bool = False
    ) -> Tuple[bool, str]:
        """
        Step 4: Predict binder-target complex using AlphaFold-Multimer
        
        Args:
            sequence_idx: Index of the designed sequence to use
            selected_models: Which AF-Multimer models to use (1-5)
            relax_prediction: Whether to relax the structure
            
        Returns:
            (success, message)
        """
        print(f"\n{'='*60}")
        print(f"STEP 4: Complex Prediction (AlphaFold-Multimer)")
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
            
            # Call AlphaFold-Multimer
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
            
            # Save ZIP file
            zip_filename = f"{self.session.project_name}_multimer_{sequence_idx+1}.zip"
            zip_path = self.output_dir / zip_filename
            with open(zip_path, 'wb') as f:
                f.write(result.content)
            self.session.complex.multimer_zip_path = str(zip_path)
            
            # Extract PDB from ZIP
            with zipfile.ZipFile(zip_path, 'r') as zf:
                names = zf.namelist()
                for name in names:
                    raw = zf.read(name)
                    text = raw.decode("utf-8", "ignore")
                    obj = json.loads(text)
                    complex_pdb = "".join(obj) if isinstance(obj, list) else str(obj)
            
            # Store results
            self.session.complex.complex_pdb = complex_pdb
            self.session.complex.docking_method = "alphafold_multimer"
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
        relax_prediction: bool = False
    ) -> Tuple[bool, str]:
        """
        Step 5: Predict complexes for multiple sequence candidates and rank them
        
        Args:
            num_candidates: Number of top sequences to evaluate
            selected_models: Which AF-Multimer models to use
            relax_prediction: Whether to relax structures
            
        Returns:
            (success, message)
        """
        print(f"\n{'='*60}")
        print(f"STEP 5: Batch Complex Prediction & Ranking")
        print(f"{'='*60}")
        
        num_candidates = min(num_candidates, len(self.session.binder.mpnn_sequences))
        rankings = []
        
        for i in range(num_candidates):
            print(f"\nEvaluating candidate {i+1}/{num_candidates}...")
            success, msg = self.run_complex_prediction(i, selected_models, relax_prediction)
            
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
        """Poll for async API results"""
        req_id = initial_response.headers.get("NVCF-REQID") or initial_response.headers.get("nvcf-reqid")
        if not req_id:
            raise Exception("No request ID in async response")
        
        print(f"⏳ Polling for results (request ID: {req_id})...")
        
        max_attempts = 180
        poll_interval = 10
        
        for attempt in range(max_attempts):
            time.sleep(poll_interval)
            
            try:
                result = requests.get(
                    f"{status_endpoint}/{req_id}",
                    headers=self.headers,
                    timeout=(10, 310)
                )
                
                if result.status_code == 200:
                    print(f"✅ Request completed after {(attempt + 1) * poll_interval}s")
                    return result.json()
                elif result.status_code == 202:
                    if (attempt + 1) % 6 == 0:  # Print every minute
                        print(f"   Still polling... ({(attempt + 1) * poll_interval}s)")
                    continue
                else:
                    raise Exception(f"Polling failed: {result.status_code} - {result.text}")
                    
            except requests.exceptions.RequestException as e:
                if attempt < max_attempts - 1:
                    continue
                raise e
        
        raise Exception(f"Timeout after {max_attempts * poll_interval // 60} minutes")
    
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
