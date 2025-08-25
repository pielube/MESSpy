#!/usr/bin/env python3
"""
MESS - Multi-Energy System Simulator Platform
A comprehensive energy system modeling and analysis platform

This wrapper provides a clean interface to the MESS simulation engine
while concealing the underlying implementation details.
"""

import os
import sys
import subprocess
import webbrowser
import time
import signal
import threading
from pathlib import Path

class MESSpyWrapper:
    """Wrapper class to launch MESS platform"""
    
    def __init__(self):
        self.app_name = "MESS - Multi-Energy System Simulator Platform"
        self.version = "2.1.0"
        self.port = 8501
        self.process = None
        self.is_running = False
        
    def print_banner(self):
        """Display the application banner"""
        banner = f"""
╔══════════════════════════════════════════════════════════════╗
║        MESS - Multi-Energy System Simulator Platform         ║
║                                                              ║
║      Multi-Energy System Simulation & Analysis               ║
║      Version: {self.version:<47}║
║                                                              ║
║   • Renewable Energy Integration                             ║
║   • Hydrogen Production & Storage                            ║
║   • Industrial End-Use Applications                          ║
║   • Economic Analysis & Optimization                         ║
║                                                              ║
╚══════════════════════════════════════════════════════════════╝
        """
        print(banner)
    
    def check_dependencies(self):
        """Check if required dependencies are available"""
        print("🔍 Checking system requirements...")
        
        # Check if streamlit_app.py exists
        if not Path("streamlit_app.py").exists():
            print("❌ Error: streamlit_app.py not found in current directory")
            return False
        
        # Check if required directories exist
        required_dirs = ["input_test_4", "techs", "data"]
        for dir_name in required_dirs:
            if not Path(dir_name).exists():
                print(f"⚠️  Warning: {dir_name} directory not found")
        
        print("✅ System requirements check completed")
        return True
    
    def start_application(self):
        """Start the MESS - Multi-Energy System Simulator Platform application"""
        print(f"\n🚀 Starting {self.app_name}...")
        print(f"📍 Application will be available at: http://localhost:{self.port}")
        print("⏳ Initializing simulation engine...")
        
        try:
            # Start the Streamlit app in a subprocess
            self.process = subprocess.Popen([
                sys.executable, "-m", "streamlit", "run", "streamlit_app.py",
                "--server.port", str(self.port),
                "--server.headless", "true",
                "--browser.gatherUsageStats", "false",
                "--server.enableCORS", "false",
                "--server.enableXsrfProtection", "false"
            ], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            
            self.is_running = True
            
            # Wait a moment for the app to start
            time.sleep(3)
            
            # Check if process is still running
            if self.process.poll() is None:
                print("✅ Application started successfully!")
                print(f"🌐 Opening browser at http://localhost:{self.port}")
                
                # Open browser
                webbrowser.open(f"http://localhost:{self.port}")
                
                print("\n📋 Application Features:")
                print("   • System Configuration & Parameter Management")
                print("   • Real-time Simulation Engine")
                print("   • Comprehensive Results Analysis")
                print("   • Advanced Data Visualization")
                print("   • Export & Reporting Capabilities")
                
                print(f"\n🔄 Application is running. Press Ctrl+C to stop.")
                
                # Keep the wrapper running
                try:
                    self.process.wait()
                except KeyboardInterrupt:
                    self.stop_application()
                    
            else:
                print("❌ Failed to start application")
                return False
                
        except Exception as e:
            print(f"❌ Error starting application: {e}")
            return False
    
    def stop_application(self):
        """Stop the MESS application"""
        if self.process and self.is_running:
            print("\n🛑 Stopping MESS application...")
            
            # Terminate the process
            self.process.terminate()
            
            # Wait for graceful shutdown
            try:
                self.process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                self.process.kill()
            
            self.is_running = False
            print("✅ Application stopped successfully")
    
    def run(self):
        """Main run method"""
        self.print_banner()
        
        if not self.check_dependencies():
            return False
        
        # Set up signal handler for graceful shutdown
        def signal_handler(signum, frame):
            self.stop_application()
            sys.exit(0)
        
        signal.signal(signal.SIGINT, signal_handler)
        signal.signal(signal.SIGTERM, signal_handler)
        
        return self.start_application()

def main():
    """Main entry point"""
    wrapper = MESSpyWrapper()
    wrapper.run()

if __name__ == "__main__":
    main()
