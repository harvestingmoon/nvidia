#!/usr/bin/env python3
"""
Quick test script to verify the NVIDIA AlphaFold2 API
"""

import requests
import json
import time

def test_alphafold2_api():
    """Test the AlphaFold2 API with a simple sequence"""
    
    # Configuration
    api_key = "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"
    endpoint = "https://api.nvcf.nvidia.com"
    function_id = "e3dfc6dd-fc27-4f0e-9ede-94412256af18"  # AlphaFold2
    
    headers = {
        "Authorization": f"Bearer {api_key}",
        "Content-Type": "application/json"
    }
    
    # Simple test sequence (insulin B-chain)
    test_sequence = "FVNQHLCGSHLVEALYLVCGERGFFYTPKT"
    
    payload = {
        "sequence": test_sequence
    }
    
    print(f"Testing AlphaFold2 API with sequence: {test_sequence}")
    print(f"Function ID: {function_id}")
    print(f"Endpoint: {endpoint}")
    
    try:
        # Make the initial request
        response = requests.post(
            f"{endpoint}/v2/nvcf/pexec/functions/{function_id}",
            headers=headers,
            json=payload,
            timeout=60
        )
        
        print(f"Initial response status: {response.status_code}")
        print(f"Response headers: {dict(response.headers)}")
        
        if response.status_code == 200:
            result = response.json()
            print("Synchronous result received:")
            print(json.dumps(result, indent=2))
            
        elif response.status_code == 202:
            result = response.json()
            print("Asynchronous processing started:")
            print(json.dumps(result, indent=2))
            
            if "reqId" in result:
                request_id = result["reqId"]
                print(f"Request ID: {request_id}")
                
                # Poll for result
                max_attempts = 30
                for attempt in range(max_attempts):
                    print(f"Polling attempt {attempt + 1}/{max_attempts}")
                    
                    poll_response = requests.get(
                        f"{endpoint}/v2/nvcf/pexec/status/{request_id}",
                        headers=headers,
                        timeout=30
                    )
                    
                    if poll_response.status_code == 200:
                        status_result = poll_response.json()
                        status = status_result.get("status", "").upper()
                        print(f"Status: {status}")
                        
                        if status == "COMPLETED":
                            # Get final result
                            result_response = requests.get(
                                f"{endpoint}/v2/nvcf/pexec/response/{request_id}",
                                headers=headers,
                                timeout=30
                            )
                            
                            if result_response.status_code == 200:
                                final_result = result_response.json()
                                print("Final result received!")
                                print("Result keys:", list(final_result.keys()))
                                
                                # Look for PDB content
                                if "pdb" in final_result:
                                    pdb_content = final_result["pdb"]
                                    print(f"PDB content length: {len(pdb_content)} characters")
                                    print("First 200 characters of PDB:")
                                    print(pdb_content[:200])
                                elif "output" in final_result:
                                    output = final_result["output"]
                                    print(f"Output type: {type(output)}")
                                    if isinstance(output, dict):
                                        print("Output keys:", list(output.keys()))
                                    elif isinstance(output, str):
                                        print(f"Output length: {len(output)} characters")
                                        print("First 200 characters:")
                                        print(output[:200])
                                else:
                                    print("Full result:")
                                    print(json.dumps(final_result, indent=2)[:1000])
                            else:
                                print(f"Failed to get result: {result_response.status_code} - {result_response.text}")
                            break
                            
                        elif status == "FAILED":
                            error_msg = status_result.get("error", "Unknown error")
                            print(f"Prediction failed: {error_msg}")
                            break
                            
                        elif status in ["PENDING", "IN_PROGRESS", "QUEUED"]:
                            print(f"Still processing... waiting 10 seconds")
                            time.sleep(10)
                            continue
                        else:
                            print(f"Unknown status: {status}")
                            break
                    else:
                        print(f"Status poll failed: {poll_response.status_code} - {poll_response.text}")
                        time.sleep(10)
                        continue
                        
                else:
                    print("Timed out waiting for result")
            else:
                print("No request ID in response")
        else:
            print(f"Request failed: {response.status_code}")
            print(f"Response: {response.text}")
            
    except requests.exceptions.RequestException as e:
        print(f"Request exception: {e}")
    except Exception as e:
        print(f"General exception: {e}")

if __name__ == "__main__":
    test_alphafold2_api()
