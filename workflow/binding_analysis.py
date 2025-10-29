"""
Binding Interface Analysis Module
Analyzes protein-protein interactions, interfaces, and binding quality
"""

import re
import numpy as np
from typing import Dict, List, Tuple, Optional, Set
from dataclasses import dataclass


@dataclass
class Atom:
    """Atom data structure"""
    serial: int
    name: str
    residue_name: str
    chain: str
    residue_number: int
    x: float
    y: float
    z: float
    element: str = ""
    
    @property
    def coords(self) -> np.ndarray:
        return np.array([self.x, self.y, self.z])


@dataclass
class Residue:
    """Residue data structure"""
    number: int
    name: str
    chain: str
    atoms: List[Atom]
    
    @property
    def ca_atom(self) -> Optional[Atom]:
        """Get CA atom of residue"""
        for atom in self.atoms:
            if atom.name.strip() == 'CA':
                return atom
        return None


def parse_pdb_content(pdb_content: str) -> List[Atom]:
    """
    Parse PDB content and extract atoms
    
    Args:
        pdb_content: PDB file content as string
        
    Returns:
        List of Atom objects
    """
    atoms = []
    
    for line in pdb_content.split('\n'):
        if line.startswith('ATOM') or line.startswith('HETATM'):
            try:
                atom = Atom(
                    serial=int(line[6:11].strip()),
                    name=line[12:16].strip(),
                    residue_name=line[17:20].strip(),
                    chain=line[21].strip() if len(line) > 21 else 'A',
                    residue_number=int(line[22:26].strip()),
                    x=float(line[30:38].strip()),
                    y=float(line[38:46].strip()),
                    z=float(line[46:54].strip()),
                    element=line[76:78].strip() if len(line) > 76 else ''
                )
                atoms.append(atom)
            except (ValueError, IndexError):
                continue
    
    return atoms


def group_atoms_by_residue(atoms: List[Atom]) -> List[Residue]:
    """Group atoms into residues"""
    residue_dict = {}
    
    for atom in atoms:
        key = (atom.chain, atom.residue_number)
        if key not in residue_dict:
            residue_dict[key] = Residue(
                number=atom.residue_number,
                name=atom.residue_name,
                chain=atom.chain,
                atoms=[]
            )
        residue_dict[key].atoms.append(atom)
    
    return list(residue_dict.values())


def calculate_distance(coord1: np.ndarray, coord2: np.ndarray) -> float:
    """Calculate Euclidean distance between two 3D coordinates"""
    return np.linalg.norm(coord1 - coord2)


def find_interface_residues(target_atoms: List[Atom], 
                           binder_atoms: List[Atom], 
                           cutoff: float = 5.0) -> Dict[str, any]:
    """
    Identify residues at the binding interface between two proteins
    
    Args:
        target_atoms: Atoms from target protein
        binder_atoms: Atoms from binder protein
        cutoff: Distance cutoff in Angstroms for interface definition
        
    Returns:
        Dictionary with interface analysis results
    """
    target_residues = group_atoms_by_residue(target_atoms)
    binder_residues = group_atoms_by_residue(binder_atoms)
    
    interface_target = set()
    interface_binder = set()
    contact_pairs = []
    all_distances = []
    
    # Find CA-CA distances below cutoff
    for target_res in target_residues:
        target_ca = target_res.ca_atom
        if not target_ca:
            continue
            
        for binder_res in binder_residues:
            binder_ca = binder_res.ca_atom
            if not binder_ca:
                continue
            
            distance = calculate_distance(target_ca.coords, binder_ca.coords)
            
            if distance < cutoff:
                interface_target.add(target_res.number)
                interface_binder.add(binder_res.number)
                contact_pairs.append({
                    'target_res': target_res.number,
                    'target_name': target_res.name,
                    'binder_res': binder_res.number,
                    'binder_name': binder_res.name,
                    'distance': distance
                })
                all_distances.append(distance)
    
    return {
        'interface_residues_target': sorted(list(interface_target)),
        'interface_residues_binder': sorted(list(interface_binder)),
        'num_contacts': len(contact_pairs),
        'contact_pairs': contact_pairs,
        'avg_distance': float(np.mean(all_distances)) if all_distances else 0.0,
        'min_distance': float(np.min(all_distances)) if all_distances else 0.0,
        'max_distance': float(np.max(all_distances)) if all_distances else 0.0
    }


def calculate_buried_surface_area(target_pdb: str, binder_pdb: str) -> float:
    """
    Estimate buried surface area (simplified calculation)
    
    Note: This is a simplified estimate. For accurate BSA, use tools like NACCESS or FreeSASA
    """
    # Placeholder - would need proper SASA calculation
    # For now, return an estimate based on interface size
    target_atoms = parse_pdb_content(target_pdb)
    binder_atoms = parse_pdb_content(binder_pdb)
    
    interface = find_interface_residues(target_atoms, binder_atoms)
    
    # Rough estimate: ~20 Å² per interface residue
    num_interface_residues = len(interface['interface_residues_target']) + len(interface['interface_residues_binder'])
    estimated_bsa = num_interface_residues * 20.0
    
    return estimated_bsa


def assess_binding_quality(interface_data: Dict) -> Dict[str, any]:
    """
    Assess the quality of a binding interface based on geometric criteria
    
    Returns:
        Dictionary with quality scores and recommendations
    """
    num_contacts = interface_data['num_contacts']
    avg_dist = interface_data['avg_distance']
    min_dist = interface_data['min_distance']
    
    quality_score = 0
    feedback = []
    warnings = []
    
    # Criterion 1: Number of contacts
    if num_contacts >= 15:
        quality_score += 40
        feedback.append("✅ Excellent number of interface contacts")
    elif num_contacts >= 10:
        quality_score += 30
        feedback.append("✅ Good number of interface contacts")
    elif num_contacts >= 5:
        quality_score += 15
        feedback.append("⚠️ Moderate number of interface contacts")
        warnings.append("Consider increasing interface size")
    else:
        feedback.append("❌ Insufficient interface contacts")
        warnings.append("Interface may be too small for stable binding")
    
    # Criterion 2: Average interface distance
    if 3.5 <= avg_dist <= 4.5:
        quality_score += 40
        feedback.append("✅ Optimal average interface distance")
    elif 3.0 <= avg_dist <= 5.0:
        quality_score += 25
        feedback.append("✅ Good interface distance")
    elif avg_dist < 3.0:
        quality_score += 10
        feedback.append("⚠️ Interface may be too tight")
        warnings.append("Check for steric clashes")
    else:
        quality_score += 15
        feedback.append("⚠️ Interface distance is suboptimal")
        warnings.append("Proteins may be too far apart")
    
    # Criterion 3: Minimum distance (clash detection)
    if min_dist >= 2.8:
        quality_score += 20
        feedback.append("✅ No severe steric clashes detected")
    elif min_dist >= 2.5:
        quality_score += 10
        feedback.append("⚠️ Possible minor steric clashes")
        warnings.append("Review interface packing")
    else:
        feedback.append("❌ Severe steric clashes detected")
        warnings.append("Significant structural optimization needed")
    
    # Get grade
    if quality_score >= 85:
        grade = "A - Excellent"
        recommendation = "This binder design looks highly promising! Proceed to experimental validation."
    elif quality_score >= 70:
        grade = "B - Good"
        recommendation = "Good binding interface. Minor optimization may improve binding."
    elif quality_score >= 50:
        grade = "C - Moderate"
        recommendation = "Interface shows potential but needs optimization. Focus on key residues."
    elif quality_score >= 30:
        grade = "D - Poor"
        recommendation = "Significant redesign recommended. Consider alternative approaches."
    else:
        grade = "F - Insufficient"
        recommendation = "Interface is inadequate. Complete redesign required."
    
    return {
        'quality_score': quality_score,
        'grade': grade,
        'feedback': feedback,
        'warnings': warnings,
        'recommendation': recommendation
    }


def combine_pdbs(target_pdb: str, binder_pdb: str, 
                 target_chain: str = 'A', binder_chain: str = 'B') -> str:
    """
    Combine two PDB structures into a single complex PDB
    
    Args:
        target_pdb: Target protein PDB content
        binder_pdb: Binder protein PDB content
        target_chain: Chain ID for target (default 'A')
        binder_chain: Chain ID for binder (default 'B')
        
    Returns:
        Combined PDB content
    """
    lines = []
    lines.append("REMARK Combined target and binder complex")
    lines.append(f"REMARK Target: Chain {target_chain}")
    lines.append(f"REMARK Binder: Chain {binder_chain}")
    
    # Process target PDB
    for line in target_pdb.split('\n'):
        if line.startswith('ATOM') or line.startswith('HETATM'):
            # Replace chain ID
            modified_line = line[:21] + target_chain + line[22:]
            lines.append(modified_line)
        elif line.startswith('TER'):
            lines.append(line)
    
    lines.append(f"TER")
    
    # Process binder PDB - renumber atoms
    atom_offset = 0
    for line in target_pdb.split('\n'):
        if line.startswith('ATOM'):
            atom_offset += 1
    
    for line in binder_pdb.split('\n'):
        if line.startswith('ATOM') or line.startswith('HETATM'):
            try:
                # Get original atom number and add offset
                orig_atom_num = int(line[6:11].strip())
                new_atom_num = orig_atom_num + atom_offset
                
                # Replace chain ID and atom number
                modified_line = (
                    line[:6] +
                    f"{new_atom_num:>5}" +
                    line[11:21] +
                    binder_chain +
                    line[22:]
                )
                lines.append(modified_line)
            except (ValueError, IndexError):
                continue
        elif line.startswith('TER'):
            lines.append(line)
    
    lines.append("END")
    
    return '\n'.join(lines)


def generate_contact_map_data(interface_data: Dict) -> Dict[str, any]:
    """
    Generate data for contact map visualization
    
    Returns:
        Dictionary with contact map data for plotting
    """
    contact_pairs = interface_data.get('contact_pairs', [])
    
    if not contact_pairs:
        return {'target_residues': [], 'binder_residues': [], 'distances': []}
    
    target_residues = [pair['target_res'] for pair in contact_pairs]
    binder_residues = [pair['binder_res'] for pair in contact_pairs]
    distances = [pair['distance'] for pair in contact_pairs]
    
    return {
        'target_residues': target_residues,
        'binder_residues': binder_residues,
        'distances': distances,
        'contact_pairs': contact_pairs
    }
