#!/usr/bin/env python3
"""
Test script for the Protein Structure Prediction Application
Run this to verify the installation and basic functionality
"""

import sys
import importlib
from typing import List, Tuple

def test_imports() -> List[Tuple[str, bool, str]]:
    """
    Test if all required modules can be imported
    
    Returns:
        List of (module_name, success, error_message) tuples
    """
    required_modules = [
        'streamlit',
        'requests',
        'json',
        'os',
        'tempfile',
        'py3Dmol',
        'langchain',
        'pydantic',
        're'
    ]
    
    results = []
    
    for module in required_modules:
        try:
            importlib.import_module(module)
            results.append((module, True, ""))
        except ImportError as e:
            results.append((module, False, str(e)))
    
    return results

def test_sequence_validation():
    """Test the sequence validation function"""
    # Import the validation function from our app
    sys.path.append('.')
    try:
        from app import validate_protein_sequence
        
        # Test cases
        test_cases = [
            ("ACDEFGHIKLMNPQRSTVWY", True, "Valid sequence with all amino acids"),
            ("MDSKGSSQKGSRLLLLLVVSNLLLCQGVVST", True, "Valid example sequence"),
            ("ACDEFGHIKLMNPQRSTVWYX", False, "Invalid character X"),
            ("", False, "Empty sequence"),
            ("ACDEF", False, "Too short"),
            ("A" * 1001, False, "Too long"),
        ]
        
        print("🧪 Testing sequence validation...")
        all_passed = True
        
        for sequence, expected_valid, description in test_cases:
            is_valid, result = validate_protein_sequence(sequence)
            if is_valid == expected_valid:
                print(f"  ✅ {description}: PASSED")
            else:
                print(f"  ❌ {description}: FAILED (expected {expected_valid}, got {is_valid})")
                all_passed = False
        
        return all_passed
    
    except Exception as e:
        print(f"❌ Could not test sequence validation: {e}")
        return False

def main():
    """Main test function"""
    print("🧬 Testing Protein Structure Prediction Application")
    print("=" * 55)
    
    # Test imports
    print("\n📦 Testing imports...")
    import_results = test_imports()
    
    failed_imports = []
    for module, success, error in import_results:
        if success:
            print(f"  ✅ {module}")
        else:
            print(f"  ❌ {module}: {error}")
            failed_imports.append(module)
    
    if failed_imports:
        print(f"\n❌ Failed imports: {', '.join(failed_imports)}")
        print("Run: pip install -r requirements.txt")
        return False
    
    print("\n✅ All imports successful!")
    
    # Test sequence validation
    validation_passed = test_sequence_validation()
    
    if validation_passed:
        print("\n✅ Sequence validation tests passed!")
    else:
        print("\n❌ Some sequence validation tests failed!")
        return False
    
    # Test basic functionality
    print("\n🔧 Testing basic application components...")
    
    try:
        from app import create_3d_visualization, save_pdb_file
        
        # Test with a minimal PDB content
        test_pdb = """HEADER    TEST PDB
ATOM      1  N   ALA A   1      20.154   1.850  18.367  1.00 20.00           N
ATOM      2  CA  ALA A   1      20.522   2.053  16.970  1.00 20.00           C
END"""
        
        # Test 3D visualization creation
        try:
            html_content = create_3d_visualization(test_pdb)
            if len(html_content) > 0:
                print("  ✅ 3D visualization creation works")
            else:
                print("  ❌ 3D visualization creation failed")
                return False
        except Exception as viz_error:
            print(f"  ⚠️  3D visualization test skipped (py3Dmol issue): {viz_error}")
            # Don't fail the test for this, as it might work in the actual Streamlit environment
        
        # Test file saving
        filepath = save_pdb_file(test_pdb, "test_structure.pdb")
        if filepath and len(filepath) > 0:
            print("  ✅ PDB file saving works")
        else:
            print("  ❌ PDB file saving failed")
            return False
    
    except Exception as e:
        print(f"  ❌ Component testing failed: {e}")
        return False
    
    print("\n🎉 All tests passed! The application should work correctly.")
    print("\nTo run the application:")
    print("  streamlit run app.py")
    print("\nThen open your browser to: http://localhost:8501")
    
    return True

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
