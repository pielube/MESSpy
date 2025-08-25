import streamlit as st
import json
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import subprocess
import os
import sys
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
from datetime import datetime
import io
import zipfile

# Set page config
st.set_page_config(
    page_title="MESS - Multi-Energy System Simulator Platform",
    page_icon="⚡",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS
st.markdown("""
<style>
    .main-header {
        font-size: 2.5rem;
        color: #1f77b4;
        text-align: center;
        margin-bottom: 2rem;
    }
    .metric-card {
        background-color: #f0f2f6;
        padding: 1rem;
        border-radius: 0.5rem;
        border-left: 4px solid #1f77b4;
    }
    .config-section {
        background-color: #ffffff;
        padding: 1.5rem;
        border-radius: 0.5rem;
        border: 1px solid #e0e0e0;
        margin-bottom: 1rem;
    }

</style>
""", unsafe_allow_html=True)

def load_config(config_file):
    """Load configuration from JSON file"""
    try:
        with open(config_file, 'r') as f:
            return json.load(f)
    except Exception as e:
        st.error(f"Error loading {config_file}: {str(e)}")
        return None

def save_config(config_data, config_file):
    """Save configuration to JSON file"""
    try:
        with open(config_file, 'w') as f:
            json.dump(config_data, f, indent=2)
        return True
    except Exception as e:
        st.error(f"Error saving {config_file}: {str(e)}")
        return False

def run_simulation():
    """Run the hybrid plant simulation"""
    try:
        result = subprocess.run([sys.executable, "run_test_4.py"], 
                              capture_output=True, text=True, cwd=os.getcwd())
        return result.returncode == 0, result.stdout, result.stderr
    except Exception as e:
        return False, "", str(e)

def display_energy_flow_diagram():
    """Display the energy flow diagram"""
    diagram_path = "input_test_4/plots/energy_flow_diagram.png"
    if os.path.exists(diagram_path):
        # Add expander for better viewing
        with st.expander("🔍 Click to expand Energy Flow Diagram", expanded=True):
            st.image(diagram_path, caption="Energy Flow Diagram", use_container_width=True)
            
            # Add download button
            with open(diagram_path, "rb") as f:
                st.download_button(
                    label="📥 Download Energy Flow Diagram",
                    data=f.read(),
                    file_name="energy_flow_diagram.png",
                    mime="image/png"
                )
    else:
        st.warning("Energy flow diagram not found. Please run the simulation first.")

def display_system_summary():
    """Display system summary metrics"""
    summary_path = "input_test_4/plots/system_summary.png"
    if os.path.exists(summary_path):
        with st.expander("🔍 Click to expand System Summary", expanded=True):
            st.image(summary_path, caption="System Summary", use_container_width=True)
            
            with open(summary_path, "rb") as f:
                st.download_button(
                    label="📥 Download System Summary",
                    data=f.read(),
                    file_name="system_summary.png",
                    mime="image/png"
                )
    else:
        st.warning("System summary not found. Please run the simulation first.")

def display_economic_analysis():
    """Display economic analysis plots"""
    analysis_path = "input_test_4/plots/economic_analysis.png"
    if os.path.exists(analysis_path):
        with st.expander("🔍 Click to expand Economic Analysis", expanded=True):
            st.image(analysis_path, caption="Economic Analysis", use_container_width=True)
            
            with open(analysis_path, "rb") as f:
                st.download_button(
                    label="📥 Download Economic Analysis",
                    data=f.read(),
                    file_name="economic_analysis.png",
                    mime="image/png"
                )
    else:
        st.warning("Economic analysis not found. Please run the simulation first.")

def display_load_profiles():
    """Display load profile plots"""
    load_path = "input_test_4/plots/load_profiles.png"
    if os.path.exists(load_path):
        with st.expander("🔍 Click to expand Load Profiles", expanded=True):
            st.image(load_path, caption="Load Profiles", use_container_width=True)
            
            with open(load_path, "rb") as f:
                st.download_button(
                    label="📥 Download Load Profiles",
                    data=f.read(),
                    file_name="load_profiles.png",
                    mime="image/png"
                )
    else:
        st.warning("Load profiles not found. Please run the simulation first.")

def display_production_profiles():
    """Display production profile plots"""
    prod_path = "input_test_4/plots/production_profiles.png"
    if os.path.exists(prod_path):
        with st.expander("🔍 Click to expand Production Profiles", expanded=True):
            st.image(prod_path, caption="Production Profiles", use_container_width=True)
            
            with open(prod_path, "rb") as f:
                st.download_button(
                    label="📥 Download Production Profiles",
                    data=f.read(),
                    file_name="production_profiles.png",
                    mime="image/png"
                )
    else:
        st.warning("Production profiles not found. Please run the simulation first.")

def display_pdf_report():
    """Display PDF report inline"""
    pdf_path = "input_test_4/plots/hybrid_plant_simulation_report.pdf"
    if os.path.exists(pdf_path):
        # File info
        file_size = os.path.getsize(pdf_path) / (1024 * 1024)  # Convert to MB
        st.info(f"📄 PDF Report Size: {file_size:.2f} MB")
        
        # Download button
        with open(pdf_path, "rb") as f:
            pdf_bytes = f.read()
        st.download_button(
            label="📥 Download PDF Report",
            data=pdf_bytes,
            file_name="hybrid_plant_simulation_report.pdf",
            mime="application/pdf"
        )
        
        st.write("**📖 PDF Report Preview:**")
        st.write("The PDF report contains comprehensive analysis of the hybrid plant simulation including:")
        st.write("- Executive summary and key findings")
        st.write("- Detailed economic analysis")
        st.write("- System performance metrics")
        st.write("- Technical specifications")
        st.write("- Recommendations and conclusions")
        
        # Display PDF inline using base64 encoding
        import base64
        pdf_base64 = base64.b64encode(pdf_bytes).decode()
        
        # Create a more responsive PDF viewer
        pdf_html = f"""
        <div style="width: 100%; height: 800px; border: 1px solid #ddd; border-radius: 5px; overflow: hidden;">
            <iframe 
                src="data:application/pdf;base64,{pdf_base64}" 
                width="100%" 
                height="100%" 
                style="border: none;"
                title="Hybrid Plant Simulation Report">
            </iframe>
        </div>
        """
        st.markdown(pdf_html, unsafe_allow_html=True)
        
        st.write("💡 **Tip**: You can zoom in/out and scroll within the PDF viewer above, or download the full report using the button above.")
        
    else:
        st.warning("PDF report not found. Please run the simulation first.")

def display_image_gallery():
    """Display all images in a gallery format"""
    plots_dir = Path("input_test_4/plots")
    if not plots_dir.exists():
        st.warning("No plots directory found. Please run the simulation first.")
        return
    
    # Get all PNG files
    image_files = list(plots_dir.glob("*.png"))
    if not image_files:
        st.warning("No image files found. Please run the simulation first.")
        return
    
    st.subheader("🖼️ Image Gallery")
    st.write("View all simulation results in one place:")
    
    # Create columns for gallery layout
    cols = st.columns(2)
    
    for i, image_file in enumerate(image_files):
        col_idx = i % 2
        with cols[col_idx]:
            st.write(f"**{image_file.stem.replace('_', ' ').title()}**")
            st.image(str(image_file), use_container_width=True)
            
            # Add download button for each image
            with open(image_file, "rb") as f:
                st.download_button(
                    label=f"📥 Download {image_file.stem}",
                    data=f.read(),
                    file_name=image_file.name,
                    mime="image/png",
                    key=f"download_{image_file.stem}"
                )
            st.divider()

def display_csv_outputs():
    """Display all CSV output files with interactive viewing"""
    outputs_dir = Path("input_test_4/outputs")
    if not outputs_dir.exists():
        st.warning("No outputs directory found. Please run the simulation first.")
        return
    
    st.write("View and download all generated CSV files:")
    
    # Create tabs for different output categories
    tab1, tab2, tab3 = st.tabs(["Intermediate Outputs", "Final Outputs", "Reconciliation Reports"])
    
    with tab1:
        st.write("**🔄 Intermediate Outputs** - Raw data processing and summary statistics")
        intermediate_dir = outputs_dir / "intermediate"
        if intermediate_dir.exists():
            csv_files = list(intermediate_dir.glob("*.csv"))
            if csv_files:
                for csv_file in csv_files:
                    with st.expander(f"📄 {csv_file.stem.replace('_', ' ').title()}", expanded=False):
                        try:
                            # Load and display CSV data
                            df = pd.read_csv(csv_file)
                            
                            # Show file info
                            col1, col2, col3 = st.columns(3)
                            with col1:
                                st.metric("Rows", len(df))
                            with col2:
                                st.metric("Columns", len(df.columns))
                            with col3:
                                file_size = csv_file.stat().st_size / 1024  # KB
                                st.metric("Size", f"{file_size:.1f} KB")
                            
                            # Display data
                            st.dataframe(df, use_container_width=True)
                            
                            # Download button
                            with open(csv_file, "rb") as f:
                                st.download_button(
                                    label=f"📥 Download {csv_file.name}",
                                    data=f.read(),
                                    file_name=csv_file.name,
                                    mime="text/csv",
                                    key=f"download_intermediate_{csv_file.stem}"
                                )
                        except Exception as e:
                            st.error(f"Error loading {csv_file.name}: {str(e)}")
            else:
                st.warning("No intermediate CSV files found.")
        else:
            st.warning("Intermediate outputs directory not found.")
    
    with tab2:
        st.write("**📈 Final Outputs** - Processed results and comprehensive analysis")
        final_dir = outputs_dir / "final"
        if final_dir.exists():
            csv_files = list(final_dir.glob("*.csv"))
            if csv_files:
                for csv_file in csv_files:
                    with st.expander(f"📄 {csv_file.stem.replace('_', ' ').title()}", expanded=False):
                        try:
                            # Load and display CSV data
                            df = pd.read_csv(csv_file)
                            
                            # Show file info
                            col1, col2, col3 = st.columns(3)
                            with col1:
                                st.metric("Rows", len(df))
                            with col2:
                                st.metric("Columns", len(df.columns))
                            with col3:
                                file_size = csv_file.stat().st_size / 1024  # KB
                                st.metric("Size", f"{file_size:.1f} KB")
                            
                            # Display data
                            st.dataframe(df, use_container_width=True)
                            
                            # Download button
                            with open(csv_file, "rb") as f:
                                st.download_button(
                                    label=f"📥 Download {csv_file.name}",
                                    data=f.read(),
                                    file_name=csv_file.name,
                                    mime="text/csv",
                                    key=f"download_final_{csv_file.stem}"
                                )
                        except Exception as e:
                            st.error(f"Error loading {csv_file.name}: {str(e)}")
            else:
                st.warning("No final CSV files found.")
        else:
            st.warning("Final outputs directory not found.")
    
    with tab3:
        st.write("**🔍 Reconciliation Reports** - Data validation and consistency checks")
        reconciliation_dir = outputs_dir / "reconciliation"
        if reconciliation_dir.exists():
            csv_files = list(reconciliation_dir.glob("*.csv"))
            if csv_files:
                for csv_file in csv_files:
                    with st.expander(f"📄 {csv_file.stem.replace('_', ' ').title()}", expanded=False):
                        try:
                            # Load and display CSV data
                            df = pd.read_csv(csv_file)
                            
                            # Show file info
                            col1, col2, col3 = st.columns(3)
                            with col1:
                                st.metric("Rows", len(df))
                            with col2:
                                st.metric("Columns", len(df.columns))
                            with col3:
                                file_size = csv_file.stat().st_size / 1024  # KB
                                st.metric("Size", f"{file_size:.1f} KB")
                            
                            # Special handling for reconciliation summary
                            if "summary_statistics" in csv_file.name:
                                st.success("✅ **Data Quality Summary**")
                                if len(df) > 0:
                                    quality_score = df.iloc[0].get('Data_Quality_Score', 0)
                                    st.metric("Overall Data Quality Score", f"{quality_score:.1f}%")
                            
                            # Display data
                            st.dataframe(df, use_container_width=True)
                            
                            # Download button
                            with open(csv_file, "rb") as f:
                                st.download_button(
                                    label=f"📥 Download {csv_file.name}",
                                    data=f.read(),
                                    file_name=csv_file.name,
                                    mime="text/csv",
                                    key=f"download_reconciliation_{csv_file.stem}"
                                )
                        except Exception as e:
                            st.error(f"Error loading {csv_file.name}: {str(e)}")
            else:
                st.warning("No reconciliation CSV files found.")
        else:
            st.warning("Reconciliation directory not found.")

def display_output_summary():
    """Display summary of all generated outputs"""
    outputs_dir = Path("input_test_4/outputs")
    if not outputs_dir.exists():
        st.warning("No outputs directory found. Please run the simulation first.")
        return
    
    
    # Count files in each directory
    summary_data = []
    
    for category in ["intermediate", "final", "reconciliation"]:
        category_dir = outputs_dir / category
        if category_dir.exists():
            csv_files = list(category_dir.glob("*.csv"))
            total_size = sum(f.stat().st_size for f in csv_files) / 1024  # KB
            
            summary_data.append({
                "Category": category.replace("_", " ").title(),
                "Files": len(csv_files),
                "Total Size (KB)": f"{total_size:.1f}",
                "Status": "✅ Available" if csv_files else "⚠️ Empty"
            })
        else:
            summary_data.append({
                "Category": category.replace("_", " ").title(),
                "Files": 0,
                "Total Size (KB)": "0.0",
                "Status": "❌ Not Found"
            })
    
    # Display summary table
    summary_df = pd.DataFrame(summary_data)
    st.dataframe(summary_df, use_container_width=True)
    
    # Show recommendations
    st.subheader("💡 Quick Actions")
    col1, col2, col3 = st.columns(3)
    
    with col1:
        if st.button("🔄 Refresh Outputs", key="refresh_outputs"):
            st.rerun()
    
    with col2:
        # Prepare a zip archive of all CSVs with folder structure preserved
        try:
            csv_files_all = list(outputs_dir.rglob("*.csv"))
            if csv_files_all:
                zip_buffer = io.BytesIO()
                with zipfile.ZipFile(zip_buffer, mode="w", compression=zipfile.ZIP_DEFLATED) as zipf:
                    for csv_path in csv_files_all:
                        arcname = csv_path.relative_to(outputs_dir)
                        zipf.write(str(csv_path), arcname=str(arcname))
                zip_buffer.seek(0)
                st.download_button(
                    label="📥 Download All CSV Files",
                    data=zip_buffer.getvalue(),
                    file_name="outputs_csv_archive.zip",
                    mime="application/zip",
                    key="download_all_csv"
                )
                st.info("Use individual download buttons in the tabs above to download specific files.")
            else:
                st.info("No CSV files found to download.")
        except Exception as e:
            st.error(f"Failed to prepare CSV archive: {e}")
    
    with col3:
        if st.button("📊 View Data Quality", key="view_quality"):
            reconciliation_dir = outputs_dir / "reconciliation"
            summary_file = reconciliation_dir / "reconciliation_summary_statistics.csv"
            if summary_file.exists():
                try:
                    df = pd.read_csv(summary_file)
                    if len(df) > 0:
                        quality_score = df.iloc[0].get('Data_Quality_Score', 0)
                        st.success(f"Data Quality Score: {quality_score:.1f}%")
                    else:
                        st.warning("No quality data available")
                except:
                    st.warning("Could not load quality data")
            else:
                st.warning("Quality report not found")

def main():
    # Header
    st.markdown('<h1 class="main-header">⚡ MESS - Multi-Energy System Simulator Platform</h1>', unsafe_allow_html=True)
    
    # Determine current page first (do not clear navigation state on reruns)
    if 'current_page' not in st.session_state or not st.session_state.current_page:
        st.session_state.current_page = "Overview"
    page = st.session_state.current_page
    
    # Sidebar
    st.sidebar.title("Navigation")
    
    # Check if outputs are available
    outputs_dir = Path("input_test_4/outputs")
    outputs_available = outputs_dir.exists() and any(outputs_dir.rglob("*.csv"))
    
    if outputs_available:
        st.sidebar.success("✅ Outputs Available")
    else:
        st.sidebar.warning("⚠️ No Outputs Found")
    
    # Navigation buttons
    st.sidebar.markdown("---")
    
    # Create navigation buttons with modern styling
    nav_buttons = [
        ("🏠 Overview", "Overview"),
        ("⚙️ Configuration", "Configuration"),
        ("📊 Results", "Results"),
        ("📚 Documentation", "Documentation")
    ]
    
    for button_text, button_page in nav_buttons:
        # Determine if this is the active page
        is_active = page == button_page
        
        # Use different button types for active vs inactive
        if is_active:
            # Active page - use primary button style
            if st.sidebar.button(
                f"**{button_text}**", 
                type="primary",
                use_container_width=True, 
                key=f"nav_{button_page.lower()}"
            ):
                st.session_state.current_page = button_page
                st.rerun()
        else:
            # Inactive page - use secondary button style
            if st.sidebar.button(
                button_text, 
                use_container_width=True, 
                key=f"nav_{button_page.lower()}"
            ):
                st.session_state.current_page = button_page
                st.rerun()
    
    if page == "Overview":
        show_overview()
    elif page == "Configuration":
        show_configuration()
    elif page == "Results":
        show_results()
    elif page == "Documentation":
        show_documentation()

def show_overview():
    """Show overview page"""
    st.header("🏭 Hybrid Plant Overview")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("""
        ### System Architecture
        The hybrid plant integrates renewable energy generation with industrial hydrogen production 
        and ammonia synthesis, optimized for Indian market conditions.
        
        **Energy Flow:**
        - Wind + Solar PV → Battery Storage
        - Battery → Electrolyzer → Hydrogen Production
        - Hydrogen → Compressor → Storage/Grid
        - Hydrogen → Ammonia Production
        """)
    
    with col2:
        st.markdown("""
        ### Key Components
        - **Renewable Energy**: Wind (50 MW) + Solar PV (40 MW)
        - **Energy Storage**: Battery (20 MWh)
        - **Hydrogen Infrastructure**: Electrolyzer (25 MW), Compressor, Storage
        - **End Use**: Ammonia Production Facility
        - **Grid Integration**: Bidirectional electricity and hydrogen exchange
        
        ### Quick Start
        1. **Configure** system parameters in Configuration tab
        2. **Run Simulation** using the button in top right
        3. **View Results** in the Results tab
        """)
    
    # Load current configuration for metrics
    studycase_config = load_config("input_test_4/studycase.json")
    if studycase_config and 'hybrid_plant' in studycase_config:
        config = studycase_config['hybrid_plant']
        
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            wind_capacity = config.get('wind', {}).get('Npower', 0) / 1000
            st.metric("Wind Capacity", f"{wind_capacity:.1f} MW")
        
        with col2:
            pv_capacity = config.get('PV', {}).get('peakP', 0) / 1000
            st.metric("Solar PV Capacity", f"{pv_capacity:.1f} MW")
        
        with col3:
            battery_capacity = config.get('battery', {}).get('capacity', 0) / 1000
            st.metric("Battery Storage", f"{battery_capacity:.1f} MWh")
        
        with col4:
            electrolyzer_capacity = config.get('electrolyzer', {}).get('Npower', 0) / 1000
            st.metric("Electrolyzer", f"{electrolyzer_capacity:.1f} MW")

def show_configuration():
    """Show configuration page"""
    # Create header with simulation button in top right
    col1, col2 = st.columns([3, 1])
    
    with col1:
        st.header("⚙️ Configuration Management")
    
    with col2:
        st.write("")  # Add some spacing
        st.write("")  # Add some spacing
        
        # Check if simulation was just completed
        simulation_completed = st.session_state.get('simulation_completed', False)
        
        if simulation_completed:
            # Show completion status instead of run button
            st.success("✅ Simulation Completed!")
            st.metric("Status", "Ready")
            
            # Show quick results summary
            st.subheader("📊 Quick Results")
            
            # Check for generated plots
            plots_dir = Path("input_test_4/plots")
            if plots_dir.exists():
                plot_files = list(plots_dir.glob("*.png"))
                if plot_files:
                    st.success(f"Generated {len(plot_files)} plot files")
                    
                    # Show key metrics
                    col1, col2 = st.columns(2)
                    with col1:
                        st.metric("Plots Generated", len(plot_files))
                    with col2:
                        st.metric("Status", "✅ Complete")
                    
                    # Add buttons for next actions
                    col1, col2 = st.columns(2)
                    with col1:
                        if st.button("📊 View Full Results", type="primary", key="view_results_after_sim"):
                            # Set session state to trigger navigation
                            st.session_state.current_page = "Results"
                            st.rerun()
                    with col2:
                        if st.button("🔄 Run Again", type="secondary", key="run_simulation_again"):
                            st.session_state.simulation_completed = False
                            st.rerun()
                else:
                    st.warning("No plot files generated")
            else:
                st.warning("No plots directory found")
        else:
            # Show run simulation button
            if st.button("🚀 Run Simulation", type="primary", use_container_width=True, key="run_simulation_config"):
                with st.spinner("Running simulation..."):
                    success, stdout, stderr = run_simulation()
                    
                    if success:
                        st.session_state.simulation_completed = True
                        st.success("✅ Simulation completed successfully!")
                        st.rerun()
                    else:
                        st.error("❌ Simulation failed!")
                        with st.expander("Error Details", expanded=True):
                            st.text(stderr)
    
    # Load configurations
    studycase_config = load_config("input_test_4/studycase.json")
    tech_costs = load_config("input_test_4/tech_cost.json")
    energy_market = load_config("input_test_4/energy_market.json")
    
    if not all([studycase_config, tech_costs, energy_market]):
        st.error("Failed to load configuration files")
        return
    
    # Show quick configuration status (only if simulation not completed)
    if not st.session_state.get('simulation_completed', False):
        st.info("💡 **Quick Tip**: Configure your system parameters below, then click 'Run Simulation' in the top right to execute the simulation.")
    else:
        st.success("🎉 **Simulation Complete!** Your results are ready. Click 'View Full Results' in the top right to see detailed analysis.")
    
    # Show current configuration summary
    if studycase_config and 'hybrid_plant' in studycase_config:
        config = studycase_config['hybrid_plant']
        
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            wind_capacity = config.get('wind', {}).get('Npower', 0) / 1000
            st.metric("Wind Power", f"{wind_capacity:.1f} MW")
        with col2:
            pv_capacity = config.get('PV', {}).get('peakP', 0) / 1000
            st.metric("Solar PV", f"{pv_capacity:.1f} MW")
        with col3:
            battery_capacity = config.get('battery', {}).get('capacity', 0) / 1000
            st.metric("Battery Storage", f"{battery_capacity:.1f} MWh")
        with col4:
            electrolyzer_capacity = config.get('electrolyzer', {}).get('Npower', 0) / 1000
            st.metric("Electrolyzer", f"{electrolyzer_capacity:.1f} MW")
    
    st.divider()
    
    # Tabs for different configuration sections
    tab1, tab2, tab3, tab4, tab5 = st.tabs(["System Components", "Technology Costs", "Energy Market", "Load Profiles", "Technical Specifications"])
    
    with tab1:
        st.subheader("System Component Configuration")
        
        if 'hybrid_plant' in studycase_config:
            config = studycase_config['hybrid_plant']
            
            # Wind Configuration
            with st.expander("🌪️ Wind Power Configuration", expanded=True):
                col1, col2, col3 = st.columns(3)
                with col1:
                    wind_power = st.number_input(
                        "Wind Power (MW)", 
                        value=float(config.get('wind', {}).get('Npower', 50000) / 1000),
                        min_value=1.0,
                        max_value=1000.0,
                        step=1.0,
                        key="wind_power_config"
                    )
                    config['wind']['Npower'] = wind_power * 1000
                
                with col2:
                    wind_cutin = st.number_input(
                        "Cut-in Speed (m/s)",
                        value=float(config.get('wind', {}).get('WScutin', 3.0)),
                        min_value=1.0,
                        max_value=10.0,
                        step=0.1,
                        key="wind_cutin_speed"
                    )
                    config['wind']['WScutin'] = wind_cutin
                
                with col3:
                    wind_rated = st.number_input(
                        "Rated Speed (m/s)",
                        value=float(config.get('wind', {}).get('WSrated', 13.0)),
                        min_value=5.0,
                        max_value=20.0,
                        step=0.5,
                        key="wind_rated_speed"
                    )
                    config['wind']['WSrated'] = wind_rated
                
                # Additional wind parameters
                col1, col2, col3 = st.columns(3)
                with col1:
                    wind_cutoff = st.number_input(
                        "Cut-off Speed (m/s)",
                        value=float(config.get('wind', {}).get('WScutoff', 25.0)),
                        min_value=15.0,
                        max_value=35.0,
                        step=0.5,
                        key="wind_cutoff_speed"
                    )
                    config['wind']['WScutoff'] = wind_cutoff
                
                with col2:
                    wind_hub_height = st.number_input(
                        "Hub Height (m)",
                        value=float(config.get('wind', {}).get('z_i', 80.0)),
                        min_value=20.0,
                        max_value=200.0,
                        step=5.0,
                        key="wind_hub_height"
                    )
                    config['wind']['z_i'] = wind_hub_height
                
                with col3:
                    wind_model = st.selectbox(
                        "Wind Model",
                        ["power_curve", "betz", "detailed"],
                        index=0,
                        key="wind_model_type"
                    )
                    config['wind']['model'] = wind_model
                
                # Wind aging and efficiency
                col1, col2, col3 = st.columns(3)
                with col1:
                    wind_efficiency = st.number_input(
                        "Overall Efficiency",
                        value=float(config.get('wind', {}).get('efficiency', 0.95)),
                        min_value=0.3,
                        max_value=1.0,
                        step=0.01,
                        key="wind_efficiency"
                    )
                    config['wind']['efficiency'] = wind_efficiency
                
                with col2:
                    wind_aging = st.checkbox(
                        "Enable Aging Effects",
                        value=config.get('wind', {}).get('ageing', False),
                        key="wind_aging_enabled"
                    )
                    config['wind']['ageing'] = wind_aging
                
                with col3:
                    wind_degradation = st.number_input(
                        "Degradation Factor (%/year)",
                        value=float(config.get('wind', {}).get('degradation factor', 1.0)),
                        min_value=0.0,
                        max_value=5.0,
                        step=0.1,
                        key="wind_degradation_factor"
                    )
                    config['wind']['degradation factor'] = wind_degradation
            
            # Solar PV Configuration
            with st.expander("☀️ Solar PV Configuration", expanded=True):
                col1, col2, col3 = st.columns(3)
                with col1:
                    pv_power = st.number_input(
                        "Solar PV Power (MW)",
                        value=float(config.get('PV', {}).get('peakP', 40000) / 1000),
                        min_value=1.0,
                        max_value=1000.0,
                        step=1.0,
                        key="pv_power_config"
                    )
                    config['PV']['peakP'] = pv_power * 1000
                
                with col2:
                    pv_tilt = st.number_input(
                        "Tilt Angle (degrees)",
                        value=float(config.get('PV', {}).get('tilt', 28)),
                        min_value=0.0,
                        max_value=90.0,
                        step=1.0,
                        key="pv_tilt_angle"
                    )
                    config['PV']['tilt'] = pv_tilt
                
                with col3:
                    pv_losses = st.number_input(
                        "System Losses (%)",
                        value=float(config.get('PV', {}).get('losses', 12)),
                        min_value=0.0,
                        max_value=30.0,
                        step=1.0,
                        key="pv_system_losses"
                    )
                    config['PV']['losses'] = pv_losses
                
                # Additional PV parameters
                col1, col2, col3 = st.columns(3)
                with col1:
                    pv_azimuth = st.number_input(
                        "Azimuth Angle (degrees)",
                        value=float(config.get('PV', {}).get('azimuth', 0)),
                        min_value=0.0,
                        max_value=360.0,
                        step=5.0,
                        key="pv_azimuth_angle"
                    )
                    config['PV']['azimuth'] = pv_azimuth
                
                with col2:
                    pv_tracking = st.selectbox(
                        "Tracking Type",
                        ["fixed", "single_axis", "dual_axis"],
                        index=0,
                        key="pv_tracking_type"
                    )
                    config['PV']['trackingtype'] = pv_tracking
                
                with col3:
                    pv_optimal_angles = st.checkbox(
                        "Auto-optimize Angles",
                        value=config.get('PV', {}).get('optimal angles', False),
                        key="pv_optimal_angles"
                    )
                    config['PV']['optimal angles'] = pv_optimal_angles
                
                # PV aging and weather data
                col1, col2, col3 = st.columns(3)
                with col1:
                    pv_aging = st.checkbox(
                        "Enable Aging Effects",
                        value=config.get('PV', {}).get('ageing', False),
                        key="pv_aging_enabled"
                    )
                    config['PV']['ageing'] = pv_aging
                
                with col2:
                    pv_degradation = st.number_input(
                        "Degradation Factor (%/year)",
                        value=float(config.get('PV', {}).get('degradation factor', 0.5)),
                        min_value=0.0,
                        max_value=2.0,
                        step=0.1,
                        key="pv_degradation_factor"
                    )
                    config['PV']['degradation factor'] = pv_degradation
                
                with col3:
                    pv_weather_source = st.selectbox(
                        "Weather Data Source",
                        ["TMY", "2023", "2022", "2021"],
                        index=0,
                        key="pv_weather_source"
                    )
                    config['PV']['serie'] = pv_weather_source
            
            # Battery Configuration
            with st.expander("🔋 Battery Storage Configuration", expanded=True):
                col1, col2, col3 = st.columns(3)
                with col1:
                    battery_capacity = st.number_input(
                        "Battery Capacity (MWh)",
                        value=float(config.get('battery', {}).get('capacity', 20000) / 1000),
                        min_value=1.0,
                        max_value=1000.0,
                        step=1.0,
                        key="battery_capacity_mwh"
                    )
                    config['battery']['capacity'] = battery_capacity * 1000
                
                with col2:
                    battery_power = st.number_input(
                        "Battery Power (MW)",
                        value=float(config.get('battery', {}).get('max power', 10000) / 1000),
                        min_value=1.0,
                        max_value=500.0,
                        step=1.0,
                        key="battery_power_mw"
                    )
                    config['battery']['max power'] = battery_power * 1000
                
                with col3:
                    battery_efficiency = st.number_input(
                        "Battery Efficiency",
                        value=float(config.get('battery', {}).get('efficiency', 0.9)),
                        min_value=0.5,
                        max_value=1.0,
                        step=0.01,
                        key="battery_efficiency"
                    )
                    config['battery']['efficiency'] = battery_efficiency
                
                # Additional battery parameters
                col1, col2, col3 = st.columns(3)
                with col1:
                    battery_charge_power = st.number_input(
                        "Max Charge Power (MW)",
                        value=float(config.get('battery', {}).get('max charging power', 10000) / 1000),
                        min_value=1.0,
                        max_value=500.0,
                        step=1.0,
                        key="battery_charge_power"
                    )
                    config['battery']['max charging power'] = battery_charge_power * 1000
                
                with col2:
                    battery_discharge_power = st.number_input(
                        "Max Discharge Power (MW)",
                        value=float(config.get('battery', {}).get('max discharging power', 10000) / 1000),
                        min_value=1.0,
                        max_value=500.0,
                        step=1.0,
                        key="battery_discharge_power"
                    )
                    config['battery']['max discharging power'] = battery_discharge_power * 1000
                
                with col3:
                    battery_dod = st.number_input(
                        "Depth of Discharge (%)",
                        value=float(config.get('battery', {}).get('depth of discharge', 90)),
                        min_value=50.0,
                        max_value=100.0,
                        step=5.0,
                        key="battery_dod"
                    )
                    config['battery']['depth of discharge'] = battery_dod
                
                # Battery aging and lifecycle
                col1, col2, col3 = st.columns(3)
                with col1:
                    battery_aging = st.checkbox(
                        "Enable Aging Effects",
                        value=config.get('battery', {}).get('ageing', False),
                        key="battery_aging_enabled"
                    )
                    config['battery']['ageing'] = battery_aging
                
                with col2:
                    battery_cycles = st.number_input(
                        "Cycle Life",
                        value=float(config.get('battery', {}).get('life cycles', 4000)),
                        min_value=1000.0,
                        max_value=10000.0,
                        step=100.0,
                        key="battery_cycle_life"
                    )
                    config['battery']['life cycles'] = battery_cycles
                
                with col3:
                    battery_end_life = st.number_input(
                        "End-of-Life Capacity (%)",
                        value=float(config.get('battery', {}).get('end life capacity', 80)),
                        min_value=50.0,
                        max_value=100.0,
                        step=5.0,
                        key="battery_end_life_capacity"
                    )
                    config['battery']['end life capacity'] = battery_end_life
            
            # Electrolyzer Configuration
            with st.expander("⚗️ Electrolyzer Configuration", expanded=True):
                col1, col2, col3 = st.columns(3)
                with col1:
                    electrolyzer_power = st.number_input(
                        "Electrolyzer Power (MW)",
                        value=float(config.get('electrolyzer', {}).get('Npower', 25000) / 1000),
                        min_value=1.0,
                        max_value=100.0,
                        step=1.0,
                        key="electrolyzer_power_mw"
                    )
                    config['electrolyzer']['Npower'] = electrolyzer_power * 1000
                
                with col2:
                    electrolyzer_modules = st.number_input(
                        "Number of Modules",
                        value=float(config.get('electrolyzer', {}).get('number of modules', 250)),
                        min_value=1.0,
                        max_value=1000.0,
                        step=10.0,
                        key="electrolyzer_modules"
                    )
                    config['electrolyzer']['number of modules'] = electrolyzer_modules
                
                with col3:
                    electrolyzer_model = st.selectbox(
                        "Stack Model",
                        ["PEM General", "SOFC", "simple"],
                        index=0,
                        key="electrolyzer_model"
                    )
                    config['electrolyzer']['stack model'] = electrolyzer_model
                
                # Additional electrolyzer parameters
                col1, col2, col3 = st.columns(3)
                with col1:
                    electrolyzer_min_load = st.number_input(
                        "Minimum Load (%)",
                        value=float(config.get('electrolyzer', {}).get('minimum_load', 10)),
                        min_value=0.0,
                        max_value=50.0,
                        step=5.0,
                        key="electrolyzer_min_load"
                    )
                    config['electrolyzer']['minimum_load'] = electrolyzer_min_load
                
                with col2:
                    electrolyzer_strategy = st.selectbox(
                        "Operation Strategy",
                        ["load_following", "baseload", "peak_shaving"],
                        index=0,
                        key="electrolyzer_strategy"
                    )
                    config['electrolyzer']['strategy'] = electrolyzer_strategy
                
                with col3:
                    electrolyzer_renewables_only = st.checkbox(
                        "Renewable-Only Operation",
                        value=config.get('electrolyzer', {}).get('only_renewables', True),
                        key="electrolyzer_renewables_only"
                    )
                    config['electrolyzer']['only_renewables'] = electrolyzer_renewables_only
                
                # Electrolyzer aging and efficiency
                col1, col2, col3 = st.columns(3)
                with col1:
                    electrolyzer_aging = st.checkbox(
                        "Enable Aging Effects",
                        value=config.get('electrolyzer', {}).get('ageing', False),
                        key="electrolyzer_aging_enabled"
                    )
                    config['electrolyzer']['ageing'] = electrolyzer_aging
                
                with col2:
                    electrolyzer_efficiency = st.number_input(
                        "System Efficiency (%)",
                        value=float(config.get('electrolyzer', {}).get('efficiency', 75)),
                        min_value=0.0,
                        max_value=90.0,
                        step=1.0,
                        key="electrolyzer_efficiency"
                    )
                    config['electrolyzer']['efficiency'] = electrolyzer_efficiency
                
                with col3:
                    electrolyzer_state = st.selectbox(
                        "System State",
                        ["on", "off", "standby"],
                        index=0,
                        key="electrolyzer_state"
                    )
                    config['electrolyzer']['state'] = electrolyzer_state
            
            # Hydrogen Storage Configuration
            with st.expander("💧 Hydrogen Storage Configuration", expanded=True):
                col1, col2, col3 = st.columns(3)
                with col1:
                    h2_storage_capacity = st.number_input(
                        "Storage Capacity (kg)",
                        value=float(config.get('hydrogen_storage', {}).get('max capacity', 5000)),
                        min_value=100.0,
                        max_value=100000.0,
                        step=100.0,
                        key="h2_storage_capacity"
                    )
                    config['hydrogen_storage']['max capacity'] = h2_storage_capacity
                
                with col2:
                    h2_storage_pressure = st.number_input(
                        "Storage Pressure (bar)",
                        value=float(config.get('hydrogen_storage', {}).get('pressure', 300)),
                        min_value=50.0,
                        max_value=1000.0,
                        step=10.0,
                        key="h2_storage_pressure"
                    )
                    config['hydrogen_storage']['pressure'] = h2_storage_pressure
                
                with col3:
                    h2_storage_owned = st.checkbox(
                        "Owned Asset",
                        value=config.get('hydrogen_storage', {}).get('owned', True),
                        key="h2_storage_owned"
                    )
                    config['hydrogen_storage']['owned'] = h2_storage_owned
                
                # Additional hydrogen storage parameters
                col1, col2, col3 = st.columns(3)
                with col1:
                    h2_storage_min_level = st.number_input(
                        "Minimum Level (%)",
                        value=float(config.get('hydrogen_storage', {}).get('min level', 10)),
                        min_value=0.0,
                        max_value=50.0,
                        step=5.0,
                        key="h2_storage_min_level"
                    )
                    config['hydrogen_storage']['min level'] = h2_storage_min_level
                
                with col2:
                    h2_storage_max_level = st.number_input(
                        "Maximum Level (%)",
                        value=float(config.get('hydrogen_storage', {}).get('max level', 95)),
                        min_value=50.0,
                        max_value=100.0,
                        step=5.0,
                        key="h2_storage_max_level"
                    )
                    config['hydrogen_storage']['max level'] = h2_storage_max_level
                
                with col3:
                    h2_storage_temperature = st.number_input(
                        "Storage Temperature (K)",
                        value=float(config.get('hydrogen_storage', {}).get('temperature', 298)),
                        min_value=200.0,
                        max_value=400.0,
                        step=5.0,
                        key="h2_storage_temperature"
                    )
                    config['hydrogen_storage']['temperature'] = h2_storage_temperature
            
            # Hydrogen Compressor Configuration
            with st.expander("🔧 Hydrogen Compressor Configuration", expanded=True):
                col1, col2, col3 = st.columns(3)
                with col1:
                    compressor_power = st.number_input(
                        "Compressor Power (kW)",
                        value=float(config.get('hydrogen_compressor', {}).get('Npower', 5000)),
                        min_value=100.0,
                        max_value=50000.0,
                        step=100.0,
                        key="compressor_power_kw"
                    )
                    config['hydrogen_compressor']['Npower'] = compressor_power
                
                with col2:
                    compressor_model = st.selectbox(
                        "Compressor Model",
                        ["simple", "normal", "with_refrigeration"],
                        index=1,
                        key="compressor_model"
                    )
                    config['hydrogen_compressor']['compressor model'] = compressor_model
                
                with col3:
                    compressor_stages = st.number_input(
                        "Number of Stages",
                        value=float(config.get('hydrogen_compressor', {}).get('n_stages', 3)),
                        min_value=1.0,
                        max_value=10.0,
                        step=1.0,
                        key="compressor_stages"
                    )
                    config['hydrogen_compressor']['n_stages'] = compressor_stages
                
                # Additional compressor parameters
                col1, col2, col3 = st.columns(3)
                with col1:
                    compressor_p_in = st.number_input(
                        "Input Pressure (bar)",
                        value=float(config.get('hydrogen_compressor', {}).get('P_in', 30)),
                        min_value=1.0,
                        max_value=100.0,
                        step=1.0,
                        key="compressor_p_in"
                    )
                    config['hydrogen_compressor']['P_in'] = compressor_p_in
                
                with col2:
                    compressor_p_out = st.number_input(
                        "Output Pressure (bar)",
                        value=float(config.get('hydrogen_compressor', {}).get('P_out', 700)),
                        min_value=100.0,
                        max_value=1000.0,
                        step=10.0,
                        key="compressor_p_out"
                    )
                    config['hydrogen_compressor']['P_out'] = compressor_p_out
                
                with col3:
                    compressor_t_in = st.number_input(
                        "Input Temperature (K)",
                        value=float(config.get('hydrogen_compressor', {}).get('T_in', 298)),
                        min_value=250.0,
                        max_value=350.0,
                        step=5.0,
                        key="compressor_t_in"
                    )
                    config['hydrogen_compressor']['T_in'] = compressor_t_in
                
                # Compressor efficiency and operation
                col1, col2, col3 = st.columns(3)
                with col1:
                    compressor_efficiency = st.number_input(
                        "Compressor Efficiency (%)",
                        value=float(config.get('hydrogen_compressor', {}).get('efficiency', 75)),
                        min_value=50.0,
                        max_value=90.0,
                        step=1.0,
                        key="compressor_efficiency"
                    )
                    config['hydrogen_compressor']['efficiency'] = compressor_efficiency
                
                with col2:
                    compressor_renewables_only = st.checkbox(
                        "Renewable-Only Operation",
                        value=config.get('hydrogen_compressor', {}).get('only_renewables', True),
                        key="compressor_renewables_only"
                    )
                    config['hydrogen_compressor']['only_renewables'] = compressor_renewables_only
                
                with col3:
                    compressor_flow_rate = st.number_input(
                        "Nominal Flow Rate (kg/h)",
                        value=float(config.get('hydrogen_compressor', {}).get('flow_rate', 100)),
                        min_value=10.0,
                        max_value=1000.0,
                        step=10.0,
                        key="compressor_flow_rate"
                    )
                    config['hydrogen_compressor']['flow_rate'] = compressor_flow_rate
            
            # Ammonia Production Configuration
            with st.expander("🏭 Ammonia Production Configuration", expanded=True):
                col1, col2, col3 = st.columns(3)
                with col1:
                    ammonia_capacity = st.number_input(
                        "Production Capacity (kg/h)",
                        value=float(config.get('ammonia_production', {}).get('capacity', 1000)),
                        min_value=100.0,
                        max_value=10000.0,
                        step=100.0,
                        key="ammonia_capacity"
                    )
                    config['ammonia_production']['capacity'] = ammonia_capacity
                
                with col2:
                    ammonia_efficiency = st.number_input(
                        "Production Efficiency (%)",
                        value=float(config.get('ammonia_production', {}).get('efficiency', 85)),
                        min_value=50.0,
                        max_value=95.0,
                        step=1.0,
                        key="ammonia_efficiency"
                    )
                    config['ammonia_production']['efficiency'] = ammonia_efficiency
                
                with col3:
                    ammonia_owned = st.checkbox(
                        "Owned Asset",
                        value=config.get('ammonia_production', {}).get('owned', True),
                        key="ammonia_owned"
                    )
                    config['ammonia_production']['owned'] = ammonia_owned
                
                # Additional ammonia production parameters
                col1, col2, col3 = st.columns(3)
                with col1:
                    ammonia_h2_consumption = st.number_input(
                        "H₂ Consumption (kg/h)",
                        value=float(config.get('ammonia_production', {}).get('h2_consumption', 177)),
                        min_value=50.0,
                        max_value=500.0,
                        step=5.0,
                        key="ammonia_h2_consumption"
                    )
                    config['ammonia_production']['h2_consumption'] = ammonia_h2_consumption
                
                with col2:
                    ammonia_n2_consumption = st.number_input(
                        "N₂ Consumption (kg/h)",
                        value=float(config.get('ammonia_production', {}).get('n2_consumption', 823)),
                        min_value=200.0,
                        max_value=2000.0,
                        step=10.0,
                        key="ammonia_n2_consumption"
                    )
                    config['ammonia_production']['n2_consumption'] = ammonia_n2_consumption
                
                with col3:
                    ammonia_heat_consumption = st.number_input(
                        "Heat Consumption (kW)",
                        value=float(config.get('ammonia_production', {}).get('heat_consumption', 500)),
                        min_value=100.0,
                        max_value=2000.0,
                        step=50.0,
                        key="ammonia_heat_consumption"
                    )
                    config['ammonia_production']['heat_consumption'] = ammonia_heat_consumption
                
                # Ammonia operation parameters
                col1, col2, col3 = st.columns(3)
                with col1:
                    ammonia_operation_hours = st.number_input(
                        "Operation Hours (h/year)",
                        value=float(config.get('ammonia_production', {}).get('operation_hours', 8000)),
                        min_value=1000.0,
                        max_value=8760.0,
                        step=100.0,
                        key="ammonia_operation_hours"
                    )
                    config['ammonia_production']['operation_hours'] = ammonia_operation_hours
                
                with col2:
                    ammonia_startup_time = st.number_input(
                        "Startup Time (hours)",
                        value=float(config.get('ammonia_production', {}).get('startup_time', 24)),
                        min_value=1.0,
                        max_value=72.0,
                        step=1.0,
                        key="ammonia_startup_time"
                    )
                    config['ammonia_production']['startup_time'] = ammonia_startup_time
                
                with col3:
                    ammonia_min_load = st.number_input(
                        "Minimum Load (%)",
                        value=float(config.get('ammonia_production', {}).get('min_load', 30)),
                        min_value=10.0,
                        max_value=80.0,
                        step=5.0,
                        key="ammonia_min_load"
                    )
                    config['ammonia_production']['min_load'] = ammonia_min_load
            
            # Heat Pump Configuration
            with st.expander("🔥 Heat Pump Configuration", expanded=False):
                col1, col2, col3 = st.columns(3)
                with col1:
                    heatpump_power = st.number_input(
                        "Heat Pump Power (kW)",
                        value=float(config.get('heat_pump', {}).get('nom Pth', 1000)),
                        min_value=100.0,
                        max_value=10000.0,
                        step=100.0,
                        key="heatpump_power"
                    )
                    config['heat_pump'] = config.get('heat_pump', {})
                    config['heat_pump']['nom Pth'] = heatpump_power
                
                with col2:
                    heatpump_type = st.selectbox(
                        "Heat Pump Type",
                        ["air_water", "water_water", "ground_water"],
                        index=0,
                        key="heatpump_type"
                    )
                    config['heat_pump']['type'] = heatpump_type
                
                with col3:
                    heatpump_cop = st.number_input(
                        "COP (Coefficient of Performance)",
                        value=float(config.get('heat_pump', {}).get('COP', 3.5)),
                        min_value=2.0,
                        max_value=6.0,
                        step=0.1,
                        key="heatpump_cop"
                    )
                    config['heat_pump']['COP'] = heatpump_cop
                
                # Additional heat pump parameters
                col1, col2, col3 = st.columns(3)
                with col1:
                    heatpump_temp_heat = st.number_input(
                        "Heating Temperature (°C)",
                        value=float(config.get('heat_pump', {}).get('t rad heat', 45)),
                        min_value=30.0,
                        max_value=80.0,
                        step=5.0,
                        key="heatpump_temp_heat"
                    )
                    config['heat_pump']['t rad heat'] = heatpump_temp_heat
                
                with col2:
                    heatpump_temp_cool = st.number_input(
                        "Cooling Temperature (°C)",
                        value=float(config.get('heat_pump', {}).get('t rad cool', 7)),
                        min_value=0.0,
                        max_value=20.0,
                        step=1.0,
                        key="heatpump_temp_cool"
                    )
                    config['heat_pump']['t rad cool'] = heatpump_temp_cool
                
                with col3:
                    heatpump_owned = st.checkbox(
                        "Owned Asset",
                        value=config.get('heat_pump', {}).get('owned', True),
                        key="heatpump_owned"
                    )
                    config['heat_pump']['owned'] = heatpump_owned
            
            # CHP Configuration
            with st.expander("⚡ CHP (Combined Heat & Power) Configuration", expanded=False):
                col1, col2, col3 = st.columns(3)
                with col1:
                    chp_electrical_power = st.number_input(
                        "Electrical Power (kW)",
                        value=float(config.get('chp', {}).get('Ppeak', 2000)),
                        min_value=100.0,
                        max_value=10000.0,
                        step=100.0,
                        key="chp_electrical_power"
                    )
                    config['chp'] = config.get('chp', {})
                    config['chp']['Ppeak'] = chp_electrical_power
                
                with col2:
                    chp_thermal_power = st.number_input(
                        "Thermal Power (kW)",
                        value=float(config.get('chp', {}).get('Qpeak', 2500)),
                        min_value=100.0,
                        max_value=15000.0,
                        step=100.0,
                        key="chp_thermal_power"
                    )
                    config['chp']['Qpeak'] = chp_thermal_power
                
                with col3:
                    chp_fuel_type = st.selectbox(
                        "Fuel Type",
                        ["natural_gas", "biogas", "hydrogen", "diesel"],
                        index=0,
                        key="chp_fuel_type"
                    )
                    config['chp']['fuel_type'] = chp_fuel_type
                
                # Additional CHP parameters
                col1, col2, col3 = st.columns(3)
                with col1:
                    chp_electrical_efficiency = st.number_input(
                        "Electrical Efficiency (%)",
                        value=float(config.get('chp', {}).get('efficiency_el', 35)),
                        min_value=20.0,
                        max_value=50.0,
                        step=1.0,
                        key="chp_electrical_efficiency"
                    )
                    config['chp']['efficiency_el'] = chp_electrical_efficiency
                
                with col2:
                    chp_thermal_efficiency = st.number_input(
                        "Thermal Efficiency (%)",
                        value=float(config.get('chp', {}).get('efficiency_th', 45)),
                        min_value=30.0,
                        max_value=60.0,
                        step=1.0,
                        key="chp_thermal_efficiency"
                    )
                    config['chp']['efficiency_th'] = chp_thermal_efficiency
                
                with col3:
                    chp_owned = st.checkbox(
                        "Owned Asset",
                        value=config.get('chp', {}).get('owned', True),
                        key="chp_owned"
                    )
                    config['chp']['owned'] = chp_owned
            
            # System-level configuration
            st.divider()
            st.subheader("🔧 System-Level Configuration")
            
            col1, col2, col3 = st.columns(3)
            with col1:
                system_simulation_hours = st.number_input(
                    "Simulation Hours",
                    value=float(config.get('system', {}).get('simulation_hours', 8760)),
                    min_value=1000.0,
                    max_value=8760.0,
                    step=100.0,
                    key="system_simulation_hours"
                )
                config['system'] = config.get('system', {})
                config['system']['simulation_hours'] = system_simulation_hours
            
            with col2:
                system_time_step = st.selectbox(
                    "Time Step",
                    ["1min", "5min", "15min", "1hour"],
                    index=0,
                    key="system_time_step"
                )
                config['system']['time_step'] = system_time_step
            
            with col3:
                system_optimization = st.checkbox(
                    "Enable Optimization",
                    value=config.get('system', {}).get('optimization', False),
                    key="system_optimization"
                )
                config['system']['optimization'] = system_optimization
            
            # Save configuration button
            st.divider()
            col1, col2, col3 = st.columns([1, 2, 1])
            with col2:
                if st.button("💾 Save System Configuration", type="primary", use_container_width=True, key="save_system_config"):
                    if save_config(studycase_config, "input_test_4/studycase.json"):
                        st.success("✅ System configuration saved successfully!")
                        st.info("💡 Configuration will be used in the next simulation run.")
                    else:
                        st.error("❌ Failed to save configuration.")
    
    with tab2:
        st.subheader("Technology Cost Configuration")
        
        # Create a form for cost updates
        with st.form("cost_config_form"):
            st.write("Update technology costs (₹/kW or ₹/kWh):")
            
            col1, col2 = st.columns(2)
            
            with col1:
                pv_cost = st.number_input(
                    "Solar PV Cost (₹/kW)",
                    value=float(tech_costs.get('PV', {}).get('cost per unit', 45000)),
                    min_value=10000.0,
                    max_value=100000.0,
                    step=1000.0
                )
                
                wind_cost = st.number_input(
                    "Wind Cost (₹/kW)",
                    value=float(tech_costs.get('wind', {}).get('cost per unit', 65000)),
                    min_value=20000.0,
                    max_value=150000.0,
                    step=1000.0
                )
                
                battery_cost = st.number_input(
                    "Battery Cost (₹/kWh)",
                    value=float(tech_costs.get('battery', {}).get('cost per unit', 75000)),
                    min_value=20000.0,
                    max_value=200000.0,
                    step=1000.0
                )
                
                electrolyzer_cost = st.number_input(
                    "Electrolyzer Cost (₹/kW)",
                    value=float(tech_costs.get('electrolyzer', {}).get('cost per unit', 135000)),
                    min_value=50000.0,
                    max_value=300000.0,
                    step=5000.0
                )
            
            with col2:
                compressor_cost = st.number_input(
                    "Hydrogen Compressor Cost (₹/kW)",
                    value=float(tech_costs.get('hydrogen_compressor', {}).get('cost per unit', 35000)),
                    min_value=10000.0,
                    max_value=100000.0,
                    step=1000.0
                )
                
                storage_cost = st.number_input(
                    "Hydrogen Storage Cost (₹/kW)",
                    value=float(tech_costs.get('hydrogen_storage', {}).get('cost per unit', 60000)),
                    min_value=10000.0,
                    max_value=150000.0,
                    step=1000.0
                )
                
                ammonia_cost = st.number_input(
                    "Ammonia Production Cost (₹/kW)",
                    value=float(tech_costs.get('ammonia_production', {}).get('cost per unit', 25000)),
                    min_value=5000.0,
                    max_value=100000.0,
                    step=1000.0
                )
            
            if st.form_submit_button("Update Costs"):
                tech_costs['PV']['cost per unit'] = pv_cost
                tech_costs['wind']['cost per unit'] = wind_cost
                tech_costs['battery']['cost per unit'] = battery_cost
                tech_costs['electrolyzer']['cost per unit'] = electrolyzer_cost
                tech_costs['hydrogen_compressor']['cost per unit'] = compressor_cost
                tech_costs['hydrogen_storage']['cost per unit'] = storage_cost
                tech_costs['ammonia_production']['cost per unit'] = ammonia_cost
                
                if save_config(tech_costs, "input_test_4/tech_cost.json"):
                    st.success("Technology costs updated successfully!")
    
    with tab3:
        st.subheader("Energy Market Configuration")
        
        with st.form("market_config_form"):
            st.write("Update energy market parameters:")
            
            col1, col2 = st.columns(2)
            
            with col1:
                electricity_purchase = st.number_input(
                    "Electricity Purchase Price (₹/kWh)",
                    value=float(energy_market.get('electricity', {}).get('purchase', 0.12)),
                    min_value=0.01,
                    max_value=1.0,
                    step=0.01
                )
                
                electricity_sale = st.number_input(
                    "Electricity Sale Price (₹/kWh)",
                    value=float(energy_market.get('electricity', {}).get('sale', 0.08)),
                    min_value=0.01,
                    max_value=1.0,
                    step=0.01
                )
                
                hydrogen_purchase = st.number_input(
                    "Hydrogen Purchase Price (₹/kg)",
                    value=float(energy_market.get('hydrogen', {}).get('purchase', 2.5)),
                    min_value=0.1,
                    max_value=10.0,
                    step=0.1
                )
                
                hydrogen_sale = st.number_input(
                    "Hydrogen Sale Price (₹/kg)",
                    value=float(energy_market.get('hydrogen', {}).get('sale', 2.0)),
                    min_value=0.1,
                    max_value=10.0,
                    step=0.1
                )
            
            with col2:
                green_h2_incentive = st.number_input(
                    "Green Hydrogen Incentive (₹/kg)",
                    value=float(energy_market.get('green_hydrogen_incentives', {}).get('value', 1.5)),
                    min_value=0.0,
                    max_value=10.0,
                    step=0.1
                )
                
                ammonia_incentive = st.number_input(
                    "Ammonia Production Incentive (₹/kg)",
                    value=float(energy_market.get('ammonia_incentives', {}).get('value', 0.05)),
                    min_value=0.0,
                    max_value=5.0,
                    step=0.01
                )
                
                interest_rate = st.number_input(
                    "Interest Rate (%)",
                    value=float(energy_market.get('interest rate', 0.085) * 100),
                    min_value=1.0,
                    max_value=20.0,
                    step=0.1
                )
            
            if st.form_submit_button("Update Market Parameters"):
                energy_market['electricity']['purchase'] = electricity_purchase
                energy_market['electricity']['sale'] = electricity_sale
                energy_market['hydrogen']['purchase'] = hydrogen_purchase
                energy_market['hydrogen']['sale'] = hydrogen_sale
                energy_market['green_hydrogen_incentives']['value'] = green_h2_incentive
                energy_market['ammonia_incentives']['value'] = ammonia_incentive
                energy_market['interest rate'] = interest_rate / 100
                
                if save_config(energy_market, "input_test_4/energy_market.json"):
                    st.success("Energy market parameters updated successfully!")
    
    with tab4:
        st.subheader("📊 Load Profile Management")
        st.info("""
        Load profiles define the energy demand patterns for the hybrid plant simulation.
        Each profile contains time-series data with 1-minute resolution for a full year.
        """)
        
        # Load profile files
        loads_dir = Path("input_test_4/loads")
        if loads_dir.exists():
            load_files = list(loads_dir.glob("*.csv"))
            if load_files:
                st.write("**📁 Available Load Profile Files:**")
                
                # Create tabs for different load types
                load_tabs = st.tabs(["Electricity Loads", "Heat Loads", "Hydrogen Loads", "File Management"])
                
                with load_tabs[0]:
                    st.write("**⚡ Electricity Load Profiles**")
                    elec_files = [f for f in load_files if 'electric' in f.name.lower()]
                    for file in elec_files:
                        with st.expander(f"📄 {file.name}", expanded=False):
                            try:
                                df = pd.read_csv(file)
                                
                                # Show file info
                                col1, col2, col3 = st.columns(3)
                                with col1:
                                    st.metric("Data Points", len(df))
                                with col2:
                                    st.metric("Columns", len(df.columns))
                                with col3:
                                    file_size = file.stat().st_size / 1024  # KB
                                    st.metric("Size", f"{file_size:.1f} KB")
                                
                                # Show data preview
                                st.write("**Data Preview (First 10 rows):**")
                                st.dataframe(df.head(10), use_container_width=True)
                                
                                # Show statistics
                                if 'electricity_demand' in df.columns:
                                    col1, col2, col3, col4 = st.columns(4)
                                    with col1:
                                        st.metric("Max Demand (kW)", f"{df['electricity_demand'].max():.1f}")
                                    with col2:
                                        st.metric("Min Demand (kW)", f"{df['electricity_demand'].min():.1f}")
                                    with col3:
                                        st.metric("Avg Demand (kW)", f"{df['electricity_demand'].mean():.1f}")
                                    with col4:
                                        st.metric("Total Energy (MWh)", f"{df['electricity_demand'].sum()/1000:.1f}")
                                
                                # Download button
                                with open(file, "rb") as f:
                                    st.download_button(
                                        label=f"📥 Download {file.name}",
                                        data=f.read(),
                                        file_name=file.name,
                                        mime="text/csv",
                                        key=f"download_elec_{file.stem}"
                                    )
                            except Exception as e:
                                st.error(f"Error loading {file.name}: {str(e)}")
                
                with load_tabs[1]:
                    st.write("**🔥 Heat Load Profiles**")
                    heat_files = [f for f in load_files if 'heat' in f.name.lower()]
                    for file in heat_files:
                        with st.expander(f"📄 {file.name}", expanded=False):
                            try:
                                df = pd.read_csv(file)
                                
                                # Show file info
                                col1, col2, col3 = st.columns(3)
                                with col1:
                                    st.metric("Data Points", len(df))
                                with col2:
                                    st.metric("Columns", len(df.columns))
                                with col3:
                                    file_size = file.stat().st_size / 1024  # KB
                                    st.metric("Size", f"{file_size:.1f} KB")
                                
                                # Show data preview
                                st.write("**Data Preview (First 10 rows):**")
                                st.dataframe(df.head(10), use_container_width=True)
                                
                                # Show statistics
                                if 'heat_demand' in df.columns:
                                    col1, col2, col3, col4 = st.columns(4)
                                    with col1:
                                        st.metric("Max Demand (kW)", f"{df['heat_demand'].max():.1f}")
                                    with col2:
                                        st.metric("Min Demand (kW)", f"{df['heat_demand'].min():.1f}")
                                    with col3:
                                        st.metric("Avg Demand (kW)", f"{df['heat_demand'].mean():.1f}")
                                    with col4:
                                        st.metric("Total Energy (MWh)", f"{df['heat_demand'].sum()/1000:.1f}")
                                
                                # Download button
                                with open(file, "rb") as f:
                                    st.download_button(
                                        label=f"📥 Download {file.name}",
                                        data=f.read(),
                                        file_name=file.name,
                                        mime="text/csv",
                                        key=f"download_heat_{file.stem}"
                                    )
                            except Exception as e:
                                st.error(f"Error loading {file.name}: {str(e)}")
                
                with load_tabs[2]:
                    st.write("**💧 Hydrogen Load Profiles**")
                    h2_files = [f for f in load_files if 'hydrogen' in f.name.lower()]
                    for file in h2_files:
                        with st.expander(f"📄 {file.name}", expanded=False):
                            try:
                                df = pd.read_csv(file)
                                
                                # Show file info
                                col1, col2, col3 = st.columns(3)
                                with col1:
                                    st.metric("Data Points", len(df))
                                with col2:
                                    st.metric("Columns", len(df.columns))
                                with col3:
                                    file_size = file.stat().st_size / 1024  # KB
                                    st.metric("Size", f"{file_size:.1f} KB")
                                
                                # Show data preview
                                st.write("**Data Preview (First 10 rows):**")
                                st.dataframe(df.head(10), use_container_width=True)
                                
                                # Show statistics
                                if 'hydrogen_demand' in df.columns:
                                    col1, col2, col3, col4 = st.columns(4)
                                    with col1:
                                        st.metric("Max Demand (kg/h)", f"{df['hydrogen_demand'].max():.1f}")
                                    with col2:
                                        st.metric("Min Demand (kg/h)", f"{df['hydrogen_demand'].min():.1f}")
                                    with col3:
                                        st.metric("Avg Demand (kg/h)", f"{df['hydrogen_demand'].mean():.1f}")
                                    with col4:
                                        st.metric("Total Demand (tons)", f"{df['hydrogen_demand'].sum()/1000:.1f}")
                                
                                # Download button
                                with open(file, "rb") as f:
                                    st.download_button(
                                        label=f"📥 Download {file.name}",
                                        data=f.read(),
                                        file_name=file.name,
                                        mime="text/csv",
                                        key=f"download_h2_{file.stem}"
                                    )
                            except Exception as e:
                                st.error(f"Error loading {file.name}: {str(e)}")
                
                with load_tabs[3]:
                    st.write("**📁 File Management**")
                    
                    # Summary table
                    file_info = []
                    for file in load_files:
                        try:
                            df = pd.read_csv(file)
                            file_size = file.stat().st_size / 1024  # KB
                            file_info.append({
                                "File": file.name,
                                "Type": "Electricity" if "electric" in file.name.lower() else 
                                       "Heat" if "heat" in file.name.lower() else 
                                       "Hydrogen" if "hydrogen" in file.name.lower() else "Other",
                                "Data Points": len(df),
                                "Columns": len(df.columns),
                                "Size (KB)": f"{file_size:.1f}",
                                "Last Modified": datetime.fromtimestamp(file.stat().st_mtime).strftime("%Y-%m-%d %H:%M")
                            })
                        except:
                            file_info.append({
                                "File": file.name,
                                "Type": "Error",
                                "Data Points": 0,
                                "Columns": 0,
                                "Size (KB)": "0.0",
                                "Last Modified": "Unknown"
                            })
                    
                    summary_df = pd.DataFrame(file_info)
                    st.dataframe(summary_df, use_container_width=True)
                    
                    # Bulk download
                    st.write("**📥 Bulk Download Options:**")
                    col1, col2 = st.columns(2)
                    
                    with col1:
                        if st.button("📥 Download All Load Files", key="download_all_loads"):
                            st.info("Use individual download buttons in the tabs above to download specific files.")
                    
                    with col2:
                        if st.button("🔄 Refresh Load Files", key="refresh_loads"):
                            st.rerun()
            else:
                st.warning("No load profile files found")
        else:
            st.warning("Load profiles directory not found")
    
    with tab5:
        st.subheader("🔬 Technical Specifications")
        st.info("""
        This section provides detailed technical specifications for all available technology components,
        based on the comprehensive technology analysis documentation.
        """)
        
        # Technology overview
        st.write("**📋 Available Technology Components**")
        
        tech_specs = {
            "Solar PV": {
                "Modeling": "Physics-based (PVGIS)",
                "Energy Carrier": "Electricity",
                "Key Parameters": ["peakP (kW)", "tilt (degrees)", "azimuth (degrees)", "losses (%)"],
                "Capabilities": ["Real weather integration", "Performance degradation", "High temporal resolution"],
                "Limitations": ["Fixed technology", "No shading effects", "Simplified degradation"]
            },
            "Wind Turbine": {
                "Modeling": "Power curve/Betz theory",
                "Energy Carrier": "Electricity",
                "Key Parameters": ["Npower (kW)", "WScutin (m/s)", "WSrated (m/s)", "WSoff (m/s)"],
                "Capabilities": ["Multiple modeling approaches", "Real wind data", "Load following"],
                "Limitations": ["Fixed technology", "No wake effects", "Simplified control"]
            },
            "Battery Storage": {
                "Modeling": "Electrochemical model",
                "Energy Carrier": "Electricity",
                "Key Parameters": ["capacity (kWh)", "max power (kW)", "efficiency", "SOC_min/max"],
                "Capabilities": ["Aging effects", "State of charge tracking", "Efficiency modeling"],
                "Limitations": ["Fixed chemistry", "Simplified aging", "No thermal effects"]
            },
            "Electrolyzer": {
                "Modeling": "Detailed PEM/Alkaline",
                "Energy Carrier": "H₂ + O₂",
                "Key Parameters": ["Npower (kW)", "number of modules", "stack model", "efficiency"],
                "Capabilities": ["Multiple technologies", "Modular design", "Chemical balance"],
                "Limitations": ["Fixed chemistry", "No dynamic response", "Simplified aging"]
            },
            "Hydrogen Storage": {
                "Modeling": "Thermodynamic model",
                "Energy Carrier": "Hydrogen",
                "Key Parameters": ["max capacity (kg)", "pressure (bar)", "owned"],
                "Capabilities": ["Pressure modeling", "Level tracking", "Capacity management"],
                "Limitations": ["Fixed technology", "No safety systems", "No leakage modeling"]
            },
            "Hydrogen Compressor": {
                "Modeling": "Thermodynamic model",
                "Energy Carrier": "Hydrogen",
                "Key Parameters": ["Npower (kW)", "compressor model", "n_stages", "P_in/P_out"],
                "Capabilities": ["Multi-stage compression", "Intercooling", "Efficiency modeling"],
                "Limitations": ["Fixed technology", "No dynamic control", "Simplified heat transfer"]
            },
            "Heat Pump": {
                "Modeling": "Detailed thermodynamic",
                "Energy Carrier": "Heat/Cooling",
                "Key Parameters": ["nom Pth (kW)", "type", "t rad heat/cool", "COP"],
                "Capabilities": ["Temperature dependence", "Thermal storage", "Dual-mode operation"],
                "Limitations": ["Limited types", "Fixed refrigerant", "No defrosting"]
            },
            "CHP Systems": {
                "Modeling": "Performance maps",
                "Energy Carrier": "Electricity + Heat",
                "Key Parameters": ["Ppeak (kW)", "Qpeak (kW)", "efficiency_el/th", "fuel_type"],
                "Capabilities": ["Performance maps", "Multiple fuels", "Load following"],
                "Limitations": ["Fixed technology", "No start-up modeling", "Simplified control"]
            }
        }
        
        # Create expandable sections for each technology
        for tech_name, specs in tech_specs.items():
            with st.expander(f"🔧 {tech_name}", expanded=False):
                col1, col2 = st.columns(2)
                
                with col1:
                    st.write("**📊 Specifications:**")
                    st.write(f"**Modeling Approach:** {specs['Modeling']}")
                    st.write(f"**Energy Carrier:** {specs['Energy Carrier']}")
                    
                    st.write("**🔑 Key Parameters:**")
                    for param in specs['Key Parameters']:
                        st.write(f"• {param}")
                
                with col2:
                    st.write("**✅ Capabilities:**")
                    for capability in specs['Capabilities']:
                        st.write(f"• {capability}")
                    
                    st.write("**⚠️ Limitations:**")
                    for limitation in specs['Limitations']:
                        st.write(f"• {limitation}")
        
        # Simulation capabilities summary
        st.divider()
        st.write("**🎯 Simulation Capabilities Summary**")
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.write("**✅ What the Simulation CAN Do:**")
            st.write("• Multi-energy systems integration")
            st.write("• High temporal resolution (sub-hourly)")
            st.write("• Real weather data integration")
            st.write("• Detailed component modeling")
            st.write("• Performance degradation")
            st.write("• Economic analysis")
            st.write("• Storage optimization")
            st.write("• Grid integration")
        
        with col2:
            st.write("**⚠️ What the Simulation CANNOT Do:**")
            st.write("• Real-time control")
            st.write("• Day-ahead optimization")
            st.write("• Grid stability analysis")
            st.write("• Component failure modeling")
            st.write("• Maintenance scheduling")
            st.write("• Environmental impact")
            st.write("• Dynamic market response")
            st.write("• Advanced control strategies")
        
        # Technical documentation link
        st.divider()
        st.write("**📚 For detailed technical specifications, refer to:**")
        st.write("• Technology Components Analysis documentation")
        st.write("• Individual component source code in `techs/` directory")
        st.write("• Component-specific parameter descriptions")



def show_results():
    """Show results page"""
    st.header("📊 Simulation Results")
    
    # Check if results exist
    plots_dir = Path("input_test_4/plots")
    if not plots_dir.exists():
        st.warning("No simulation results found. Please run the simulation first.")
        return
    
    # Results tabs
    tab1, tab2, tab3, tab4, tab5, tab6, tab7, tab8, tab9 = st.tabs([
        "Output Summary", "CSV Output Files", "Energy Flow", "System Summary", "Economic Analysis", "Load Profiles", "Production Profiles", "PDF Report", "Image Gallery"
    ])
    
    with tab1:
        st.subheader("📊 Output Summary")
        display_output_summary()
    
    with tab2:
        st.subheader("📊 CSV Output Files")
        display_csv_outputs()
    
    with tab3:
        st.subheader("Energy Flow Diagram")
        display_energy_flow_diagram()
    
    with tab4:
        st.subheader("System Summary")
        display_system_summary()
    
    with tab5:
        st.subheader("Economic Analysis")
        display_economic_analysis()
    
    with tab6:
        st.subheader("Load Profiles")
        display_load_profiles()
    
    with tab7:
        st.subheader("Production Profiles")
        display_production_profiles()
    
    with tab8:
        st.subheader("PDF Report")
        display_pdf_report()
    
    with tab9:
        st.subheader("Image Gallery")
        display_image_gallery()
    
    # Download section
    st.subheader("📥 Download All Results")
    
    if plots_dir.exists():
        # Show file information
        st.write("**📊 Available Files:**")
        
        # Get all files
        all_files = list(plots_dir.glob("*"))
        if all_files:
            file_info = []
            for file in all_files:
                size_mb = file.stat().st_size / (1024 * 1024)
                file_info.append({
                    "File": file.name,
                    "Type": file.suffix.upper(),
                    "Size (MB)": f"{size_mb:.2f}"
                })
            
            # Display file table
            import pandas as pd
            df = pd.DataFrame(file_info)
            st.dataframe(df, use_container_width=True)
            
            # Download buttons
            st.write("**📥 Download Options:**")
            col1, col2 = st.columns(2)
            
            with col1:
                st.write("**Individual Files:**")
                for file in all_files:
                    with open(file, "rb") as f:
                        st.download_button(
                            label=f"📥 {file.name}",
                            data=f.read(),
                            file_name=file.name,
                            mime="application/octet-stream",
                            key=f"download_{file.name}"
                        )
            
            with col2:
                st.write("**Quick Actions:**")
                if st.button("🔄 Refresh Results", key="refresh_results"):
                    st.rerun()
                
                st.write("💡 **Tip**: Use the individual tabs above for detailed viewing of each result type.")
        else:
            st.warning("No result files found. Please run the simulation first.")

def show_documentation():
    """Show documentation page"""
    st.header("📚 Documentation")
    
    # Create tabs for different documentation sections
    tab1, tab2, tab3, tab4 = st.tabs(["Model Inputs", "Technology Components", "Output Generation", "Quick Guide"])
    
    with tab1:
        st.subheader("📋 Model Inputs Documentation")
        st.markdown("""
        This section provides a comprehensive overview of all inputs being taken by the model, 
        including what variables are considered, what is assumed implicitly, and what is not used.
        """)
        
        # Load and display the Model Inputs Documentation
        try:
            with open("input_test_4/MODEL_INPUTS_DOCUMENTATION.md", "r", encoding="utf-8") as f:
                model_docs = f.read()
            
            # Convert markdown to HTML for better display
            try:
                import markdown
                html_content = markdown.markdown(model_docs, extensions=['tables', 'fenced_code'])
            except ImportError:
                # Fallback if markdown library is not available
                st.warning("Markdown library not available. Displaying raw markdown.")
                st.markdown(model_docs)
                return
            
            # Display with custom CSS for better formatting
            st.markdown("""
            <style>
            .documentation-content {
                background-color: #f8f9fa;
                padding: 20px;
                border-radius: 10px;
                border-left: 4px solid #1f77b4;
                margin: 10px 0;
            }
            .documentation-content table {
                border-collapse: collapse;
                width: 100%;
                margin: 10px 0;
            }
            .documentation-content th, .documentation-content td {
                border: 1px solid #ddd;
                padding: 8px;
                text-align: left;
            }
            .documentation-content th {
                background-color: #1f77b4;
                color: white;
            }
            .documentation-content h1, .documentation-content h2, .documentation-content h3 {
                color: #1f77b4;
                margin-top: 20px;
                margin-bottom: 10px;
            }
            .documentation-content code {
                background-color: #e9ecef;
                padding: 2px 4px;
                border-radius: 3px;
                font-family: 'Courier New', monospace;
            }
            </style>
            """, unsafe_allow_html=True)
            
            st.markdown(f'<div class="documentation-content">{html_content}</div>', unsafe_allow_html=True)
            
            # Add download button for the documentation
            with open("input_test_4/MODEL_INPUTS_DOCUMENTATION.md", "r", encoding="utf-8") as f:
                st.download_button(
                    label="📥 Download Model Inputs Documentation (PDF)",
                    data=f.read(),
                    file_name="MODEL_INPUTS_DOCUMENTATION.md",
                    mime="text/markdown"
                )
                
        except FileNotFoundError:
            st.error("Model Inputs Documentation file not found. Please ensure the file exists in the input_test_4 directory.")
        except Exception as e:
            st.error(f"Error loading Model Inputs Documentation: {str(e)}")
    
    with tab2:
        st.subheader("🔧 Technology Components Analysis")
        st.markdown("""
        This section provides a comprehensive analysis of all technology components, 
        their inputs/outputs, capabilities, and limitations to give a clear picture 
        of what the simulation can and cannot do.
        """)
        
        # Load and display the Technology Components Analysis
        try:
            with open("input_test_4/TECHNOLOGY_COMPONENTS_ANALYSIS.md", "r", encoding="utf-8") as f:
                tech_docs = f.read()
            
            # Convert markdown to HTML for better display
            try:
                import markdown
                html_content = markdown.markdown(tech_docs, extensions=['tables', 'fenced_code'])
            except ImportError:
                # Fallback if markdown library is not available
                st.warning("Markdown library not available. Displaying raw markdown.")
                st.markdown(tech_docs)
                return
            
            # Display with custom CSS for better formatting
            st.markdown("""
            <style>
            .tech-documentation-content {
                background-color: #f8f9fa;
                padding: 20px;
                border-radius: 10px;
                border-left: 4px solid #28a745;
                margin: 10px 0;
            }
            .tech-documentation-content table {
                border-collapse: collapse;
                width: 100%;
                margin: 10px 0;
            }
            .tech-documentation-content th, .tech-documentation-content td {
                border: 1px solid #ddd;
                padding: 8px;
                text-align: left;
            }
            .tech-documentation-content th {
                background-color: #28a745;
                color: white;
            }
            .tech-documentation-content h1, .tech-documentation-content h2, .tech-documentation-content h3 {
                color: #28a745;
                margin-top: 20px;
                margin-bottom: 10px;
            }
            .tech-documentation-content code {
                background-color: #e9ecef;
                padding: 2px 4px;
                border-radius: 3px;
                font-family: 'Courier New', monospace;
            }
            </style>
            """, unsafe_allow_html=True)
            
            st.markdown(f'<div class="tech-documentation-content">{html_content}</div>', unsafe_allow_html=True)
            
            # Add download button for the documentation
            with open("input_test_4/TECHNOLOGY_COMPONENTS_ANALYSIS.md", "r", encoding="utf-8") as f:
                st.download_button(
                    label="📥 Download Technology Components Analysis (PDF)",
                    data=f.read(),
                    file_name="TECHNOLOGY_COMPONENTS_ANALYSIS.md",
                    mime="text/markdown"
                )
                
        except FileNotFoundError:
            st.error("Technology Components Analysis file not found. Please ensure the file exists in the input_test_4 directory.")
        except Exception as e:
            st.error(f"Error loading Technology Components Analysis: {str(e)}")
    
    with tab3:
        st.subheader("📊 Output Generation System")
        st.markdown("""
        This section provides comprehensive documentation of the output generation system 
        that creates intermediate and final CSV files, along with a detailed reconciliation 
        exercise to validate data consistency and robustness.
        """)
        
        # Load and display the Output Generation Documentation
        try:
            with open("input_test_4/OUTPUT_GENERATION_DOCUMENTATION.md", "r", encoding="utf-8") as f:
                output_docs = f.read()
            
            # Convert markdown to HTML for better display
            try:
                import markdown
                html_content = markdown.markdown(output_docs, extensions=['tables', 'fenced_code'])
            except ImportError:
                # Fallback if markdown library is not available
                st.warning("Markdown library not available. Displaying raw markdown.")
                st.markdown(output_docs)
                return
            
            # Display with custom CSS for better formatting
            st.markdown("""
            <style>
            .output-documentation-content {
                background-color: #f8f9fa;
                padding: 20px;
                border-radius: 10px;
                border-left: 4px solid #ffc107;
                margin: 10px 0;
            }
            .output-documentation-content table {
                border-collapse: collapse;
                width: 100%;
                margin: 10px 0;
            }
            .output-documentation-content th, .output-documentation-content td {
                border: 1px solid #ddd;
                padding: 8px;
                text-align: left;
            }
            .output-documentation-content th {
                background-color: #ffc107;
                color: white;
            }
            .output-documentation-content h1, .output-documentation-content h2, .output-documentation-content h3 {
                color: #ffc107;
                margin-top: 20px;
                margin-bottom: 10px;
            }
            .output-documentation-content code {
                background-color: #e9ecef;
                padding: 2px 4px;
                border-radius: 3px;
                font-family: 'Courier New', monospace;
            }
            </style>
            """, unsafe_allow_html=True)
            
            st.markdown(f'<div class="output-documentation-content">{html_content}</div>', unsafe_allow_html=True)
            
            # Add download button for the documentation
            with open("input_test_4/OUTPUT_GENERATION_DOCUMENTATION.md", "r", encoding="utf-8") as f:
                st.download_button(
                    label="📥 Download Output Generation Documentation (PDF)",
                    data=f.read(),
                    file_name="OUTPUT_GENERATION_DOCUMENTATION.md",
                    mime="text/markdown"
                )
                
        except FileNotFoundError:
            st.error("Output Generation Documentation file not found. Please ensure the file exists in the input_test_4 directory.")
        except Exception as e:
            st.error(f"Error loading Output Generation Documentation: {str(e)}")
    
    with tab4:
        st.subheader("🚀 Quick Start Guide")
        st.markdown("""
        ### Hybrid Plant Simulation Quick Guide
        
        This application provides an interactive interface for configuring and running 
        a hybrid renewable energy and hydrogen production plant simulation.
        
        #### System Components
        
        **🌪️ Wind Power**
        - Capacity: Configurable wind power generation
        - Technology: Power curve model with configurable parameters
        - Location: Optimized for Indian wind conditions
        
        **☀️ Solar PV**
        - Capacity: Configurable solar power generation
        - Technology: Fixed tilt with tracking options
        - Location: Optimized for Indian solar conditions
        
        **🔋 Battery Storage**
        - Capacity: Configurable energy storage
        - Technology: Lithium-ion battery with efficiency modeling
        - Purpose: Energy time-shifting and grid stability
        
        **⚗️ Electrolyzer**
        - Capacity: Configurable hydrogen production
        - Technology: PEM electrolyzer with modular design
        - Purpose: Green hydrogen production from renewable electricity
        
        **💧 Hydrogen Infrastructure**
        - Compressor: Multi-stage hydrogen compression
        - Storage: High-pressure hydrogen storage
        - Purpose: Hydrogen handling and distribution
        
        **🏭 Ammonia Production**
        - Capacity: Industrial ammonia synthesis facility
        - Technology: Green ammonia production from hydrogen
        - Purpose: End-use application and energy carrier
        
        #### Configuration Files
        
        - **studycase.json**: System component configuration
        - **tech_cost.json**: Technology cost parameters
        - **energy_market.json**: Energy prices and incentives
        - **Load profiles**: CSV files with demand patterns
        
        #### Simulation Outputs
        
        - Energy flow diagrams
        - Economic analysis
        - System performance metrics
        - Load and production profiles
        - PDF summary report
        
        #### Usage Instructions
        
        1. **Configuration**: Modify system parameters in the Configuration tab
        2. **Simulation**: Run the simulation in the Simulation tab
        3. **Results**: View and download results in the Results tab
        4. **Documentation**: Reference the detailed documentation tabs above
        
        #### Key Features
        
        ✅ **127 Explicit Input Parameters** across 6 major categories  
        ✅ **16 Technology Components** with detailed modeling  
        ✅ **Real Weather Data Integration** via PVGIS  
        ✅ **Comprehensive Economic Analysis** with Indian market data  
        ✅ **High Temporal Resolution** simulation capability  
        ✅ **Multi-Energy System** modeling (electricity, heat, hydrogen, oxygen)  
        
        #### Model Capabilities
        
        **What the simulation CAN do:**
        - System design and component sizing
        - Economic feasibility analysis
        - Technology comparison studies
        - Long-term strategic planning
        - Policy impact assessment
        - Multi-year performance analysis
        
        **What the simulation CANNOT do:**
        - Real-time control and operation
        - Day-ahead optimization
        - Grid stability analysis
        - Component failure modeling
        - Maintenance scheduling
        
        #### Getting Started
        
        1. **Review Documentation**: Start with the detailed documentation tabs
        2. **Configure System**: Use the Configuration tab to set up your system
        3. **Run Simulation**: Execute the simulation in the Simulation tab
        4. **Analyze Results**: Review comprehensive results in the Results tab
        5. **Download Reports**: Save all analysis and documentation for further use
        """)
        
        # Add a summary metrics section
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Input Parameters", "127", "Explicit parameters")
        with col2:
            st.metric("Technology Components", "16", "Available technologies")
        with col3:
            st.metric("Energy Carriers", "4", "Electricity, Heat, H₂, O₂")
        
        # Add a capabilities summary
        st.subheader("📊 Model Summary")
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("""
            **✅ Strengths:**
            - Comprehensive multi-energy modeling
            - High technical detail
            - Real weather integration
            - Detailed economic analysis
            - Long-term analysis capability
            - Modular component design
            """)
        
        with col2:
            st.markdown("""
            **⚠️ Limitations:**
            - No real-time control
            - Simplified operations
            - Fixed technology options
            - No environmental impact
            - Basic grid integration
            """)



if __name__ == "__main__":
    main()
