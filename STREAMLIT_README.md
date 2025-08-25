# ⚡ Hybrid Plant Simulation Dashboard

A comprehensive Streamlit web application for configuring, running, and analyzing hybrid renewable energy and hydrogen production plant simulations.

## 🚀 Features

### 📊 Interactive Dashboard
- **Overview**: System architecture and key metrics
- **Configuration**: Modify system parameters, costs, and market conditions
- **Simulation**: Run the hybrid plant simulation with real-time feedback
- **Results**: View and download simulation outputs
- **Documentation**: Comprehensive system documentation

### ⚙️ Configuration Management
- **System Components**: Wind, Solar PV, Battery, Electrolyzer, Hydrogen Infrastructure
- **Technology Costs**: Update capital and O&M costs for all components
- **Energy Market**: Modify electricity prices, hydrogen prices, and incentives
- **Load Profiles**: Information about demand patterns

### 📈 Results Visualization
- Energy flow diagrams
- Economic analysis plots
- System summary metrics
- Load and production profiles
- Downloadable results and PDF reports

## 🛠️ Installation

1. **Install Dependencies**:
   ```bash
   pip install -r streamlit_requirements.txt
   ```

2. **Ensure Simulation Files**:
   - Make sure `run_test_4.py` is in the current directory
   - Ensure all configuration files exist in `input_test_4/`

## 🎯 Usage

### Starting the Application
```bash
streamlit run streamlit_app.py
```

The application will open in your default web browser at `http://localhost:8501`

### Navigation

#### 1. **Overview Page**
- View system architecture and energy flow
- See current system metrics and capacities
- Understand the hybrid plant concept

#### 2. **Configuration Page**
- **System Components Tab**: Modify wind, solar, battery, electrolyzer, and hydrogen storage parameters
- **Technology Costs Tab**: Update capital and O&M costs for all technologies
- **Energy Market Tab**: Adjust electricity prices, hydrogen prices, and incentives
- **Load Profiles Tab**: Information about demand patterns

#### 3. **Simulation Page**
- Check configuration file status
- View current system parameters
- Run the simulation with one click
- Monitor simulation progress and results

#### 4. **Results Page**
- View generated plots and diagrams
- Download simulation results
- Access PDF reports

#### 5. **Documentation Page**
- Comprehensive system documentation
- Component descriptions and parameters
- Usage instructions and best practices

## 📋 Configuration Parameters

### System Components
- **Wind Power**: Capacity, cut-in/rated speeds, efficiency
- **Solar PV**: Capacity, tilt angle, system losses
- **Battery Storage**: Capacity, power rating, efficiency
- **Electrolyzer**: Power rating, number of modules, stack model
- **Hydrogen Storage**: Capacity, pressure

### Technology Costs
- Capital costs (₹/kW or ₹/kWh)
- O&M costs (₹/kW/year)
- Replacement rates and lifetimes

### Energy Market
- Electricity purchase/sale prices
- Hydrogen purchase/sale prices
- Green hydrogen incentives
- Ammonia production incentives
- Interest rates

## 🔧 Technical Details

### File Structure
```
├── streamlit_app.py              # Main Streamlit application
├── streamlit_requirements.txt    # Python dependencies
├── run_test_4.py                # Simulation script
├── input_test_4/                # Configuration directory
│   ├── studycase.json           # System configuration
│   ├── tech_cost.json           # Technology costs
│   ├── energy_market.json       # Market parameters
│   ├── loads/                   # Load profile files
│   └── plots/                   # Generated plots
└── STREAMLIT_README.md          # This file
```

### Dependencies
- **Streamlit**: Web application framework
- **Pandas**: Data manipulation
- **Plotly**: Interactive visualizations
- **Matplotlib**: Plot generation
- **NumPy**: Numerical computations

## 🎨 Customization

### Adding New Components
1. Update the configuration forms in `show_configuration()`
2. Add corresponding parameters to JSON files
3. Update the simulation script if needed

### Modifying Visualizations
1. Edit the display functions in the results section
2. Add new plot types as needed
3. Customize CSS styles in the main function

### Extending Functionality
1. Add new pages to the sidebar navigation
2. Create new configuration sections
3. Implement additional analysis features

## 🐛 Troubleshooting

### Common Issues

1. **Configuration Files Not Found**
   - Ensure all JSON files exist in `input_test_4/`
   - Check file permissions

2. **Simulation Fails**
   - Verify Python dependencies are installed
   - Check simulation script (`run_test_4.py`) exists
   - Review error messages in the simulation output

3. **Plots Not Displaying**
   - Run the simulation first to generate plots
   - Check that `input_test_4/plots/` directory exists
   - Verify plot files are generated successfully

### Getting Help
- Check the Documentation page in the application
- Review simulation output for error messages
- Ensure all dependencies are properly installed

## 📄 License

This application is part of the MESSpy hybrid plant simulation framework.

## 🤝 Contributing

To contribute to this application:
1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Test thoroughly
5. Submit a pull request

---

**Note**: This application requires the underlying MESSpy simulation framework to be properly configured and functional.
