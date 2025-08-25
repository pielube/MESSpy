HYBRID PLANT SIMULATION - PLOTS DIRECTORY
==========================================

This directory contains all plots generated during the hybrid plant simulation.
All plots are saved in high-resolution PNG format (300 DPI) for later review and analysis.

GENERATED PLOTS:
================

1. LOAD_PROFILES.PNG
   - 2x2 subplot layout showing 24-hour load profiles
   - Electricity demand (with steel/ammonia breakdown)
   - Heat demand (with steel/ammonia breakdown)  
   - Hydrogen demand (with steel/ammonia breakdown)
   - Combined load comparison (all energy types)
   - File size: ~2-3 MB

2. PRODUCTION_PROFILES.PNG
   - 2x1 subplot layout showing renewable generation
   - Wind power generation profile
   - Solar PV generation profile
   - 24-hour production patterns
   - File size: ~1-2 MB

3. ENERGY_FLOW_DIAGRAM.PNG
   - Visual representation of energy flows in the hybrid plant
   - Component layout with color coding
   - Energy flow arrows showing system connectivity
   - Wind → Battery → Electrolyzer → H2 Compressor → Grid/H2 Storage → Ammonia Production
   - Solar → Battery → Electrolyzer → H2 Compressor → Grid/H2 Storage → Ammonia Production
   - File size: ~1-2 MB

4. ECONOMIC_ANALYSIS.PNG
   - 2x2 subplot layout showing economic breakdown
   - Capital costs by technology
   - O&M costs by technology
   - Energy prices comparison
   - Incentives and subsidies breakdown
   - File size: ~2-3 MB

5. SYSTEM_SUMMARY.PNG
   - 1x2 subplot layout showing system overview
   - System component capacities (bar chart)
   - Capital cost distribution (pie chart)
   - Technology cost breakdown
   - File size: ~1-2 MB

PLOT FEATURES:
==============
- High resolution: 300 DPI for print quality
- Professional styling with grid lines and legends
- Color-coded components and energy flows
- Value labels on charts for easy reading
- Consistent formatting across all plots
- Optimized for both screen viewing and printing

USAGE:
======
- Run the simulation: python run_test_4.py
- All plots are automatically generated and saved
- Plots can be opened in any image viewer
- Suitable for reports, presentations, and analysis
- Can be imported into documents or presentations

TECHNICAL DETAILS:
==================
- Format: PNG (Portable Network Graphics)
- Resolution: 300 DPI (dots per inch)
- Color space: RGB
- Compression: Lossless
- File naming: Descriptive and consistent
- Directory structure: Organized by plot type

PLOT GENERATION TIMING:
=======================
- Load profiles: ~2-3 seconds
- Production profiles: ~1-2 seconds  
- Energy flow diagram: ~1-2 seconds
- Economic analysis: ~2-3 seconds
- System summary: ~1-2 seconds
- Total generation time: ~7-12 seconds

These plots provide comprehensive visualization of the hybrid plant's:
- Energy demand patterns and industrial breakdown
- Renewable energy generation profiles
- System architecture and energy flows
- Economic feasibility and cost structure
- Overall system performance and capacity
