#!/usr/bin/env python3
"""
Output Generator for MESSpy Simulation
Generates intermediate and final output files in CSV format for analysis and reconciliation
"""

import os
import pandas as pd
import numpy as np
from datetime import datetime, timedelta
import json

class OutputGenerator:
    def __init__(self, output_dir="input_test_4/outputs"):
        """Initialize the output generator"""
        self.output_dir = output_dir
        self.intermediate_dir = os.path.join(output_dir, "intermediate")
        self.final_dir = os.path.join(output_dir, "final")
        self.reconciliation_dir = os.path.join(output_dir, "reconciliation")
        
        # Create output directories
        os.makedirs(self.output_dir, exist_ok=True)
        os.makedirs(self.intermediate_dir, exist_ok=True)
        os.makedirs(self.final_dir, exist_ok=True)
        os.makedirs(self.reconciliation_dir, exist_ok=True)
        
        # Initialize data storage
        self.simulation_data = {}
        self.intermediate_results = {}
        self.final_results = {}
        
    def generate_timestamp_series(self, start_date="2023-01-01", periods=8760, freq="1H"):
        """Generate timestamp series for the simulation period"""
        return pd.date_range(start=start_date, periods=periods, freq=freq)
    
    def save_intermediate_outputs(self, loads_data, production_data, studycase_config):
        """Generate and save intermediate output files"""
        print("📊 Generating intermediate output files...")
        
        # Generate timestamp series
        timestamps = self.generate_timestamp_series()
        
        # 1. Load Profiles Summary
        self._generate_load_summary(loads_data, timestamps)
        
        # 2. Production Data Summary
        self._generate_production_summary(production_data, timestamps)
        
        # 3. System Configuration Summary
        self._generate_system_config_summary(studycase_config)
        
        # 4. Weather Data Summary (if available)
        self._generate_weather_summary(timestamps)
        
        print(f"✅ Intermediate outputs saved to: {self.intermediate_dir}")
    
    def _generate_load_summary(self, loads_data, timestamps):
        """Generate load profiles summary"""
        load_summary = pd.DataFrame(index=timestamps)
        
        for load_file, data in loads_data.items():
            if 'electricity_demand' in data.columns:
                load_summary[f'{load_file}_electricity_kW'] = data['electricity_demand']
            if 'heat_demand' in data.columns:
                load_summary[f'{load_file}_heat_kW'] = data['heat_demand']
            if 'hydrogen_demand' in data.columns:
                load_summary[f'{load_file}_hydrogen_kg_h'] = data['hydrogen_demand']
        
        # Add summary statistics
        load_summary['total_electricity_kW'] = load_summary.filter(like='electricity').sum(axis=1)
        load_summary['total_heat_kW'] = load_summary.filter(like='heat').sum(axis=1)
        load_summary['total_hydrogen_kg_h'] = load_summary.filter(like='hydrogen').sum(axis=1)
        
        # Save to CSV
        output_file = os.path.join(self.intermediate_dir, "load_profiles_summary.csv")
        load_summary.to_csv(output_file)
        
        # Generate load statistics
        load_stats = pd.DataFrame({
            'Metric': ['Max Demand', 'Min Demand', 'Average Demand', 'Total Energy'],
            'Electricity (kW)': [
                load_summary['total_electricity_kW'].max(),
                load_summary['total_electricity_kW'].min(),
                load_summary['total_electricity_kW'].mean(),
                load_summary['total_electricity_kW'].sum()
            ],
            'Heat (kW)': [
                load_summary['total_heat_kW'].max(),
                load_summary['total_heat_kW'].min(),
                load_summary['total_heat_kW'].mean(),
                load_summary['total_heat_kW'].sum()
            ],
            'Hydrogen (kg/h)': [
                load_summary['total_hydrogen_kg_h'].max(),
                load_summary['total_hydrogen_kg_h'].min(),
                load_summary['total_hydrogen_kg_h'].mean(),
                load_summary['total_hydrogen_kg_h'].sum()
            ]
        })
        
        stats_file = os.path.join(self.intermediate_dir, "load_statistics.csv")
        load_stats.to_csv(stats_file, index=False)
    
    def _generate_production_summary(self, production_data, timestamps):
        """Generate production data summary"""
        if not production_data:
            return
            
        production_summary = pd.DataFrame(index=timestamps)
        
        # Process solar data
        solar_sources = [key for key in production_data.keys() if 'SS' in key]
        if solar_sources:
            solar_data = pd.DataFrame(index=timestamps)
            for source in solar_sources[:4]:  # Limit to first 4 sources
                data = production_data[source]
                if len(data) >= len(timestamps):
                    solar_data[f'{source}_solar_power_A_m2'] = data['solar_power'].head(len(timestamps)).values
                else:
                    # Pad with zeros if data is shorter
                    padded_data = np.pad(data['solar_power'].values, 
                                       (0, len(timestamps) - len(data)), 
                                       mode='constant')
                    solar_data[f'{source}_solar_power_A_m2'] = padded_data
            
            production_summary = pd.concat([production_summary, solar_data], axis=1)
        
        # Process wind data
        wind_sources = [key for key in production_data.keys() if 'WS' in key]
        if wind_sources:
            wind_data = pd.DataFrame(index=timestamps)
            for source in wind_sources[:4]:  # Limit to first 4 sources
                data = production_data[source]
                if len(data) >= len(timestamps):
                    wind_data[f'{source}_wind_power_kW'] = data['wind_power'].head(len(timestamps)).values
                else:
                    # Pad with zeros if data is shorter
                    padded_data = np.pad(data['wind_power'].values, 
                                       (0, len(timestamps) - len(data)), 
                                       mode='constant')
                    wind_data[f'{source}_wind_power_kW'] = padded_data
            
            production_summary = pd.concat([production_summary, wind_data], axis=1)
        
        # Save to CSV
        output_file = os.path.join(self.intermediate_dir, "production_data_summary.csv")
        production_summary.to_csv(output_file)
        
        # Generate production statistics
        if not production_summary.empty:
            prod_stats = []
            for col in production_summary.columns:
                if 'solar' in col:
                    energy_type = 'Solar (A/m²)'
                elif 'wind' in col:
                    energy_type = 'Wind (kW)'
                else:
                    energy_type = 'Other'
                
                prod_stats.append({
                    'Source': col,
                    'Type': energy_type,
                    'Max': production_summary[col].max(),
                    'Min': production_summary[col].min(),
                    'Average': production_summary[col].mean(),
                    'Total': production_summary[col].sum()
                })
            
            prod_stats_df = pd.DataFrame(prod_stats)
            stats_file = os.path.join(self.intermediate_dir, "production_statistics.csv")
            prod_stats_df.to_csv(stats_file, index=False)
    
    def _generate_system_config_summary(self, studycase_config):
        """Generate system configuration summary"""
        if not studycase_config or 'hybrid_plant' not in studycase_config:
            return
            
        config_data = []
        hybrid_config = studycase_config['hybrid_plant']
        
        for component, config in hybrid_config.items():
            if isinstance(config, dict):
                # Extract capacity information
                capacity = None
                capacity_key = None
                
                for key in ['capacity', 'Npower', 'peakP', 'max capacity', 'max power']:
                    if key in config and isinstance(config[key], (int, float)):
                        capacity = config[key]
                        capacity_key = key
                        break
                
                config_data.append({
                    'Component': component,
                    'Capacity_Value': capacity,
                    'Capacity_Unit': capacity_key,
                    'Priority': config.get('priority', 'N/A'),
                    'Owned': config.get('owned', 'N/A'),
                    'Model': config.get('model', 'N/A'),
                    'Strategy': config.get('strategy', 'N/A')
                })
        
        config_df = pd.DataFrame(config_data)
        output_file = os.path.join(self.intermediate_dir, "system_configuration.csv")
        config_df.to_csv(output_file, index=False)
    
    def _generate_weather_summary(self, timestamps):
        """Generate weather data summary if available"""
        weather_file = "input_test_4/weather/TMY_general.csv"
        if os.path.exists(weather_file):
            try:
                weather_data = pd.read_csv(weather_file)
                weather_summary = pd.DataFrame(index=timestamps)
                
                # Add weather data columns
                if len(weather_data) >= len(timestamps):
                    weather_summary['temperature_C'] = weather_data['temp_air'].head(len(timestamps)).values
                    weather_summary['relative_humidity_percent'] = weather_data['relative_humidity'].head(len(timestamps)).values
                    weather_summary['global_horizontal_irradiance_W_m2'] = weather_data['ghi'].head(len(timestamps)).values
                    weather_summary['wind_speed_m_s'] = weather_data['wind_speed'].head(len(timestamps)).values
                    weather_summary['wind_direction_degrees'] = weather_data['wind_direction'].head(len(timestamps)).values
                    weather_summary['pressure_Pa'] = weather_data['pressure'].head(len(timestamps)).values
                
                output_file = os.path.join(self.intermediate_dir, "weather_data_summary.csv")
                weather_summary.to_csv(output_file)
                
                # Generate weather statistics
                weather_stats = pd.DataFrame({
                    'Metric': ['Max', 'Min', 'Average', 'Standard Deviation'],
                    'Temperature (°C)': [
                        weather_summary['temperature_C'].max(),
                        weather_summary['temperature_C'].min(),
                        weather_summary['temperature_C'].mean(),
                        weather_summary['temperature_C'].std()
                    ],
                    'Wind Speed (m/s)': [
                        weather_summary['wind_speed_m_s'].max(),
                        weather_summary['wind_speed_m_s'].min(),
                        weather_summary['wind_speed_m_s'].mean(),
                        weather_summary['wind_speed_m_s'].std()
                    ],
                    'Solar Irradiance (W/m²)': [
                        weather_summary['global_horizontal_irradiance_W_m2'].max(),
                        weather_summary['global_horizontal_irradiance_W_m2'].min(),
                        weather_summary['global_horizontal_irradiance_W_m2'].mean(),
                        weather_summary['global_horizontal_irradiance_W_m2'].std()
                    ]
                })
                
                stats_file = os.path.join(self.intermediate_dir, "weather_statistics.csv")
                weather_stats.to_csv(stats_file, index=False)
                
            except Exception as e:
                print(f"Warning: Could not process weather data: {e}")
    
    def generate_final_outputs(self, studycase_config, tech_costs, energy_market, loads_data, production_data):
        """Generate and save final output files"""
        print("📊 Generating final output files...")
        
        # Generate timestamp series
        timestamps = self.generate_timestamp_series()
        
        # 1. Economic Analysis Summary
        self._generate_economic_summary(tech_costs, energy_market)
        
        # 2. System Performance Summary
        self._generate_performance_summary(studycase_config, loads_data, production_data, timestamps)
        
        # 3. Energy Balance Summary
        self._generate_energy_balance_summary(loads_data, production_data, timestamps)
        
        # 4. Component Analysis Summary
        self._generate_component_analysis_summary(studycase_config, tech_costs)
        
        print(f"✅ Final outputs saved to: {self.final_dir}")
    
    def _generate_economic_summary(self, tech_costs, energy_market):
        """Generate economic analysis summary"""
        # Technology costs summary
        tech_cost_data = []
        for tech, costs in tech_costs.items():
            tech_cost_data.append({
                'Technology': tech,
                'Cost_per_Unit': costs.get('cost per unit', 0),
                'OeM_Cost': costs.get('OeM', 0),
                'Replacement_Rate': costs.get('replacement', {}).get('rate', 0),
                'Replacement_Years': costs.get('replacement', {}).get('years', 0)
            })
        
        tech_cost_df = pd.DataFrame(tech_cost_data)
        output_file = os.path.join(self.final_dir, "technology_costs_summary.csv")
        tech_cost_df.to_csv(output_file, index=False)
        
        # Energy market summary
        market_data = []
        for energy_type, prices in energy_market.items():
            if isinstance(prices, dict):
                market_data.append({
                    'Energy_Type': energy_type,
                    'Purchase_Price': prices.get('purchase', 0),
                    'Sale_Price': prices.get('sale', 0),
                    'Incentive_Value': prices.get('value', 0) if 'incentive' in energy_type else 0
                })
        
        market_df = pd.DataFrame(market_data)
        output_file = os.path.join(self.final_dir, "energy_market_summary.csv")
        market_df.to_csv(output_file, index=False)
    
    def _generate_performance_summary(self, studycase_config, loads_data, production_data, timestamps):
        """Generate system performance summary"""
        performance_data = pd.DataFrame(index=timestamps)
        
        # Add load data
        if loads_data:
            if 'electric_load_hybrid.csv' in loads_data:
                performance_data['electricity_demand_kW'] = loads_data['electric_load_hybrid.csv']['electricity_demand']
            if 'heat_load_hybrid.csv' in loads_data:
                performance_data['heat_demand_kW'] = loads_data['heat_load_hybrid.csv']['heat_demand']
            if 'hydrogen_demand_hybrid.csv' in loads_data:
                performance_data['hydrogen_demand_kg_h'] = loads_data['hydrogen_demand_hybrid.csv']['hydrogen_demand']
        
        # Add production data (simplified)
        if production_data:
            solar_sources = [key for key in production_data.keys() if 'SS' in key]
            wind_sources = [key for key in production_data.keys() if 'WS' in key]
            
            if solar_sources:
                # Use first solar source as representative
                solar_data = production_data[solar_sources[0]]
                if len(solar_data) >= len(timestamps):
                    performance_data['solar_production_A_m2'] = solar_data['solar_power'].head(len(timestamps)).values
            
            if wind_sources:
                # Use first wind source as representative
                wind_data = production_data[wind_sources[0]]
                if len(wind_data) >= len(timestamps):
                    performance_data['wind_production_kW'] = wind_data['wind_power'].head(len(timestamps)).values
        
        # Save performance data
        output_file = os.path.join(self.final_dir, "system_performance_summary.csv")
        performance_data.to_csv(output_file)
        
        # Generate performance statistics
        if not performance_data.empty:
            perf_stats = []
            for col in performance_data.columns:
                perf_stats.append({
                    'Metric': col,
                    'Max': performance_data[col].max(),
                    'Min': performance_data[col].min(),
                    'Average': performance_data[col].mean(),
                    'Total': performance_data[col].sum(),
                    'Standard_Deviation': performance_data[col].std()
                })
            
            perf_stats_df = pd.DataFrame(perf_stats)
            stats_file = os.path.join(self.final_dir, "performance_statistics.csv")
            perf_stats_df.to_csv(stats_file, index=False)
    
    def _generate_energy_balance_summary(self, loads_data, production_data, timestamps):
        """Generate energy balance summary"""
        balance_data = pd.DataFrame(index=timestamps)
        
        # Calculate energy balance
        if loads_data and production_data:
            # Electricity balance
            if 'electric_load_hybrid.csv' in loads_data:
                balance_data['electricity_demand_kW'] = loads_data['electric_load_hybrid.csv']['electricity_demand']
            
            # Estimate production (simplified)
            solar_sources = [key for key in production_data.keys() if 'SS' in key]
            if solar_sources:
                solar_data = production_data[solar_sources[0]]
                if len(solar_data) >= len(timestamps):
                    # Convert solar irradiance to power (simplified)
                    solar_power = solar_data['solar_power'].head(len(timestamps)).values * 0.1  # Simplified conversion
                    balance_data['solar_power_kW'] = solar_power
            
            # Calculate net electricity
            if 'electricity_demand_kW' in balance_data.columns and 'solar_power_kW' in balance_data.columns:
                balance_data['net_electricity_kW'] = balance_data['solar_power_kW'] - balance_data['electricity_demand_kW']
                balance_data['electricity_surplus_kW'] = balance_data['net_electricity_kW'].clip(lower=0)
                balance_data['electricity_deficit_kW'] = balance_data['net_electricity_kW'].clip(upper=0).abs()
        
        # Save energy balance
        output_file = os.path.join(self.final_dir, "energy_balance_summary.csv")
        balance_data.to_csv(output_file)
    
    def _generate_component_analysis_summary(self, studycase_config, tech_costs):
        """Generate component analysis summary"""
        if not studycase_config or 'hybrid_plant' not in studycase_config:
            return
            
        component_data = []
        hybrid_config = studycase_config['hybrid_plant']
        
        for component, config in hybrid_config.items():
            if isinstance(config, dict):
                # Extract capacity
                capacity = None
                for key in ['capacity', 'Npower', 'peakP', 'max capacity', 'max power']:
                    if key in config and isinstance(config[key], (int, float)):
                        capacity = config[key]
                        break
                
                # Get cost information
                cost_per_unit = tech_costs.get(component, {}).get('cost per unit', 0)
                total_cost = capacity * cost_per_unit if capacity and cost_per_unit else 0
                
                component_data.append({
                    'Component': component,
                    'Capacity': capacity,
                    'Cost_per_Unit': cost_per_unit,
                    'Total_Cost': total_cost,
                    'Priority': config.get('priority', 'N/A'),
                    'Owned': config.get('owned', 'N/A'),
                    'Model': config.get('model', 'N/A')
                })
        
        component_df = pd.DataFrame(component_data)
        output_file = os.path.join(self.final_dir, "component_analysis_summary.csv")
        component_df.to_csv(output_file, index=False)
    
    def generate_reconciliation_report(self, loads_data, production_data, studycase_config, tech_costs, energy_market):
        """Generate reconciliation report to verify data consistency"""
        print("🔍 Generating reconciliation report...")
        
        reconciliation_data = []
        
        # 1. Data completeness check
        reconciliation_data.append({
            'Check_Type': 'Data_Completeness',
            'Load_Profiles': len(loads_data),
            'Production_Sources': len(production_data),
            'System_Components': len(studycase_config.get('hybrid_plant', {})) if studycase_config else 0,
            'Technology_Costs': len(tech_costs),
            'Energy_Market_Parameters': len(energy_market)
        })
        
        # 2. Capacity consistency check
        total_capacity = 0
        if studycase_config and 'hybrid_plant' in studycase_config:
            hybrid_config = studycase_config['hybrid_plant']
            for component, config in hybrid_config.items():
                if isinstance(config, dict):
                    capacity = None
                    for key in ['capacity', 'Npower', 'peakP', 'max capacity', 'max power']:
                        if key in config and isinstance(config[key], (int, float)):
                            capacity = config[key]
                            break
                    if capacity:
                        total_capacity += capacity
        
        reconciliation_data.append({
            'Check_Type': 'Capacity_Consistency',
            'Total_System_Capacity_kW': total_capacity,
            'Renewable_Capacity_kW': total_capacity,  # Simplified
            'Storage_Capacity_kWh': 0,  # To be calculated
            'Conversion_Capacity_kW': 0  # To be calculated
        })
        
        # 3. Load profile consistency
        max_electricity_demand = 0
        max_heat_demand = 0
        max_hydrogen_demand = 0
        
        if loads_data:
            if 'electric_load_hybrid.csv' in loads_data:
                max_electricity_demand = loads_data['electric_load_hybrid.csv']['electricity_demand'].max()
            if 'heat_load_hybrid.csv' in loads_data:
                max_heat_demand = loads_data['heat_load_hybrid.csv']['heat_demand'].max()
            if 'hydrogen_demand_hybrid.csv' in loads_data:
                max_hydrogen_demand = loads_data['hydrogen_demand_hybrid.csv']['hydrogen_demand'].max()
        
        reconciliation_data.append({
            'Check_Type': 'Load_Profile_Consistency',
            'Max_Electricity_Demand_kW': max_electricity_demand,
            'Max_Heat_Demand_kW': max_heat_demand,
            'Max_Hydrogen_Demand_kg_h': max_hydrogen_demand,
            'Total_Annual_Electricity_MWh': max_electricity_demand * 8760 / 1000 if max_electricity_demand else 0
        })
        
        # 4. Economic consistency
        total_capital_cost = 0
        total_om_cost = 0
        
        for tech, costs in tech_costs.items():
            if 'cost per unit' in costs and isinstance(costs['cost per unit'], (int, float)):
                total_capital_cost += costs['cost per unit']
            if 'OeM' in costs and isinstance(costs['OeM'], (int, float)):
                total_om_cost += costs['OeM']
        
        reconciliation_data.append({
            'Check_Type': 'Economic_Consistency',
            'Total_Capital_Cost': total_capital_cost,
            'Total_OM_Cost': total_om_cost,
            'Average_Electricity_Price': energy_market.get('electricity', {}).get('purchase', 0),
            'Average_Hydrogen_Price': energy_market.get('hydrogen', {}).get('purchase', 0)
        })
        
        # Save reconciliation report
        reconciliation_df = pd.DataFrame(reconciliation_data)
        output_file = os.path.join(self.reconciliation_dir, "reconciliation_report.csv")
        reconciliation_df.to_csv(output_file, index=False)
        
        # Generate summary statistics
        summary_stats = {
            'Total_Load_Profiles': len(loads_data),
            'Total_Production_Sources': len(production_data),
            'Total_System_Components': len(studycase_config.get('hybrid_plant', {})) if studycase_config else 0,
            'Total_Technology_Costs': len(tech_costs),
            'Total_Energy_Market_Parameters': len(energy_market),
            'Max_Electricity_Demand_kW': max_electricity_demand,
            'Max_Heat_Demand_kW': max_heat_demand,
            'Max_Hydrogen_Demand_kg_h': max_hydrogen_demand,
            'Total_System_Capacity_kW': total_capacity,
            'Total_Capital_Cost': total_capital_cost
        }
        
        summary_df = pd.DataFrame([summary_stats])
        summary_file = os.path.join(self.reconciliation_dir, "reconciliation_summary.csv")
        summary_df.to_csv(summary_file, index=False)
        
        print(f"✅ Reconciliation report saved to: {self.reconciliation_dir}")
        
        return reconciliation_df
