# MESSpy - Hybrid Plant Simulation Platform

## Overview

MESSpy is a comprehensive energy system modeling and analysis platform designed for hybrid plant simulation, renewable energy integration, and industrial hydrogen production analysis. The platform provides advanced capabilities for system design, economic analysis, and performance optimization.

## Features

### 🏭 System Components
- **Renewable Energy Integration**: Wind and solar power generation modeling
- **Energy Storage Systems**: Battery and hydrogen storage optimization
- **Hydrogen Infrastructure**: Electrolyzer, compressor, and storage systems
- **Industrial Applications**: Ammonia production and other end-use facilities
- **Grid Integration**: Bidirectional energy exchange capabilities

### 📊 Analysis Capabilities
- **Real-time Simulation**: High-resolution temporal analysis
- **Economic Modeling**: Comprehensive cost-benefit analysis
- **Performance Optimization**: System efficiency and capacity planning
- **Data Visualization**: Advanced charts and interactive dashboards
- **Export & Reporting**: PDF reports and CSV data export

### 🔧 Technical Specifications
- **127 Explicit Parameters**: Comprehensive system configuration
- **16 Technology Components**: Detailed component modeling
- **Multi-Energy Systems**: Electricity, heat, hydrogen, and oxygen
- **High Temporal Resolution**: Sub-hourly simulation capability
- **Real Weather Integration**: PVGIS and meteorological data

## Quick Start

### Windows Users
1. Double-click `launch_messpy.bat`
2. The application will start automatically
3. Your browser will open to the MESSpy interface

### Unix/Linux/macOS Users
1. Open terminal in the MESSpy directory
2. Run: `./launch_messpy.sh`
3. The application will start automatically
4. Your browser will open to the MESSpy interface

### Manual Launch
1. Open terminal/command prompt in the MESSpy directory
2. Run: `python wrapper.py`
3. The application will start automatically

## Application Structure

### 📋 Overview Page
- System architecture and component overview
- Quick start guide and navigation
- Current system capacity metrics

### ⚙️ Configuration Page
- System component parameter management
- Technology cost configuration
- Energy market settings
- Load profile management
- Technical specifications

### 📊 Results Page
- Simulation output analysis
- Interactive data visualization
- Export capabilities
- Comprehensive reporting

### 📚 Documentation Page
- Model inputs documentation
- Technology component analysis
- Output generation documentation
- User guide and tutorials

## System Requirements

### Minimum Requirements
- **Operating System**: Windows 10+, macOS 10.14+, or Linux
- **Python**: Version 3.8 or higher
- **Memory**: 4 GB RAM
- **Storage**: 2 GB available disk space

### Recommended Requirements
- **Operating System**: Windows 11, macOS 12+, or Ubuntu 20.04+
- **Python**: Version 3.9 or higher
- **Memory**: 8 GB RAM or higher
- **Storage**: 5 GB available disk space
- **Browser**: Chrome, Firefox, Safari, or Edge (latest version)

## Installation

### Automatic Installation
The MESSpy platform includes all necessary dependencies and will install them automatically on first run.

### Manual Installation (if needed)
```bash
pip install -r requirements.txt
```

## Usage Guide

### 1. System Configuration
- Navigate to the Configuration page
- Set up your hybrid plant components
- Configure technology costs and market parameters
- Review technical specifications

### 2. Simulation Execution
- Click the "Run Simulation" button in the top right
- Monitor simulation progress
- Review quick results summary

### 3. Results Analysis
- Navigate to the Results page
- Explore different analysis tabs
- Download reports and data files
- Generate custom visualizations

### 4. Documentation
- Access comprehensive documentation
- Review model inputs and assumptions
- Understand technology capabilities
- Reference user guides

## File Structure

```
MESSpy/
├── wrapper.py                 # Main application launcher
├── launch_messpy.bat         # Windows launch script
├── launch_messpy.sh          # Unix/Linux/macOS launch script
├── streamlit_app.py          # Main application interface
├── run_test_4.py            # Simulation engine
├── input_test_4/            # Configuration and data files
│   ├── studycase.json       # System configuration
│   ├── tech_cost.json       # Technology costs
│   ├── energy_market.json   # Market parameters
│   ├── loads/               # Load profile data
│   ├── outputs/             # Simulation outputs
│   └── plots/               # Generated visualizations
├── techs/                   # Technology component models
├── data/                    # Meteorological and production data
└── requirements.txt         # Python dependencies
```

## Support

### Documentation
- Comprehensive documentation is available within the application
- Technical specifications and model details are provided
- User guides and tutorials are included

### Troubleshooting
1. **Application won't start**: Check Python installation and dependencies
2. **Browser doesn't open**: Manually navigate to http://localhost:8501
3. **Simulation errors**: Review configuration files and data inputs
4. **Performance issues**: Check system requirements and available memory

## Version Information

- **Current Version**: 2.1.0
- **Release Date**: 2024
- **Compatibility**: Python 3.8+, Modern Web Browsers

## License

This software is proprietary and confidential. Unauthorized copying, distribution, or use is strictly prohibited.

---

**MESSpy - Advanced Energy System Modeling Platform**
*Empowering the future of renewable energy integration*
