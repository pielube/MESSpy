#!/usr/bin/env python3
"""
Test script for the Streamlit application
"""

import json
import os
from pathlib import Path

def test_configuration_files():
    """Test if all required configuration files exist"""
    print("Testing configuration files...")
    
    required_files = [
        "input_test_4/studycase.json",
        "input_test_4/tech_cost.json", 
        "input_test_4/energy_market.json"
    ]
    
    for file_path in required_files:
        if os.path.exists(file_path):
            print(f"✅ {file_path} - Found")
            # Test if JSON is valid
            try:
                with open(file_path, 'r') as f:
                    json.load(f)
                print(f"   ✅ Valid JSON")
            except json.JSONDecodeError as e:
                print(f"   ❌ Invalid JSON: {e}")
        else:
            print(f"❌ {file_path} - Missing")

def test_load_profiles():
    """Test if load profile files exist"""
    print("\nTesting load profile files...")
    
    loads_dir = Path("input_test_4/loads")
    if loads_dir.exists():
        load_files = list(loads_dir.glob("*.csv"))
        if load_files:
            print(f"✅ Found {len(load_files)} load profile files:")
            for file in load_files:
                print(f"   - {file.name}")
        else:
            print("❌ No load profile files found")
    else:
        print("❌ Loads directory not found")

def test_simulation_script():
    """Test if simulation script exists"""
    print("\nTesting simulation script...")
    
    if os.path.exists("run_test_4.py"):
        print("✅ run_test_4.py - Found")
    else:
        print("❌ run_test_4.py - Missing")

def test_streamlit_app():
    """Test if Streamlit app exists"""
    print("\nTesting Streamlit application...")
    
    if os.path.exists("streamlit_app.py"):
        print("✅ streamlit_app.py - Found")
    else:
        print("❌ streamlit_app.py - Missing")

def main():
    print("🧪 Testing Streamlit Application Setup")
    print("=" * 50)
    
    test_configuration_files()
    test_load_profiles()
    test_simulation_script()
    test_streamlit_app()
    
    print("\n" + "=" * 50)
    print("🎯 Test Summary")
    print("If all tests pass, you can run:")
    print("   streamlit run streamlit_app.py")
    print("\nThe application will be available at:")
    print("   http://localhost:8501")

if __name__ == "__main__":
    main()
