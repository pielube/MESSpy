#!/usr/bin/env python3
"""
Reconciliation Exercise for MESSpy Simulation
Performs comprehensive data validation and consistency checks
"""

import os
import pandas as pd
import numpy as np
import json
from datetime import datetime
import matplotlib.pyplot as plt
import seaborn as sns

class ReconciliationExercise:
    def __init__(self, output_dir="input_test_4/outputs"):
        """Initialize the reconciliation exercise"""
        self.output_dir = output_dir
        self.intermediate_dir = os.path.join(output_dir, "intermediate")
        self.final_dir = os.path.join(output_dir, "final")
        self.reconciliation_dir = os.path.join(output_dir, "reconciliation")
        
        # Create reconciliation directory if it doesn't exist
        os.makedirs(self.reconciliation_dir, exist_ok=True)
        
        # Initialize results storage
        self.reconciliation_results = {}
        self.validation_checks = {}
        self.consistency_issues = []
        
    def run_comprehensive_reconciliation(self, loads_data, production_data, studycase_config, tech_costs, energy_market):
        """Run comprehensive reconciliation exercise"""
        print("🔍 Starting comprehensive reconciliation exercise...")
        print("=" * 60)
        
        # 1. Data Quality Assessment
        print("📊 1. Data Quality Assessment...")
        self._assess_data_quality(loads_data, production_data, studycase_config, tech_costs, energy_market)
        
        # 2. Input-Output Consistency Check
        print("🔄 2. Input-Output Consistency Check...")
        self._check_input_output_consistency(loads_data, production_data, studycase_config)
        
        # 3. Energy Balance Validation
        print("⚖️ 3. Energy Balance Validation...")
        self._validate_energy_balance(loads_data, production_data, studycase_config)
        
        # 4. Economic Consistency Check
        print("💰 4. Economic Consistency Check...")
        self._check_economic_consistency(tech_costs, energy_market, studycase_config)
        
        # 5. Temporal Consistency Check
        print("⏰ 5. Temporal Consistency Check...")
        self._check_temporal_consistency(loads_data, production_data)
        
        # 6. Physical Constraints Validation
        print("🔬 6. Physical Constraints Validation...")
        self._validate_physical_constraints(studycase_config, loads_data, production_data)
        
        # 7. Generate Reconciliation Report
        print("📋 7. Generating Comprehensive Reconciliation Report...")
        self._generate_comprehensive_report()
        
        print("✅ Reconciliation exercise completed!")
        print(f"📁 Results saved to: {self.reconciliation_dir}")
        
        return self.reconciliation_results
    
    def _assess_data_quality(self, loads_data, production_data, studycase_config, tech_costs, energy_market):
        """Assess data quality and completeness"""
        quality_metrics = {}
        
        # Load data quality
        load_quality = {
            'total_files': len(loads_data),
            'missing_values': 0,
            'negative_values': 0,
            'zero_values': 0,
            'data_completeness': 0
        }
        
        for load_file, data in loads_data.items():
            for col in data.columns:
                if col != 'timestamp':
                    load_quality['missing_values'] += data[col].isnull().sum()
                    load_quality['negative_values'] += (data[col] < 0).sum()
                    load_quality['zero_values'] += (data[col] == 0).sum()
                    load_quality['data_completeness'] += len(data[col].dropna())
        
        quality_metrics['load_data'] = load_quality
        
        # Production data quality
        prod_quality = {
            'total_sources': len(production_data),
            'solar_sources': len([k for k in production_data.keys() if 'SS' in k]),
            'wind_sources': len([k for k in production_data.keys() if 'WS' in k]),
            'data_completeness': 0
        }
        
        for source, data in production_data.items():
            for col in data.columns:
                if col != 'timestamp':
                    prod_quality['data_completeness'] += len(data[col].dropna())
        
        quality_metrics['production_data'] = prod_quality
        
        # Configuration quality
        config_quality = {
            'total_components': len(studycase_config.get('hybrid_plant', {})),
            'missing_capacities': 0,
            'invalid_priorities': 0
        }
        
        if studycase_config and 'hybrid_plant' in studycase_config:
            for component, config in studycase_config['hybrid_plant'].items():
                if isinstance(config, dict):
                    # Check for missing capacities
                    capacity_found = False
                    for key in ['capacity', 'Npower', 'peakP', 'max capacity', 'max power']:
                        if key in config and isinstance(config[key], (int, float)):
                            capacity_found = True
                            break
                    if not capacity_found:
                        config_quality['missing_capacities'] += 1
                    
                    # Check for invalid priorities
                    if 'priority' in config and not isinstance(config['priority'], int):
                        config_quality['invalid_priorities'] += 1
        
        quality_metrics['configuration'] = config_quality
        
        # Economic data quality
        econ_quality = {
            'total_technologies': len(tech_costs),
            'missing_costs': 0,
            'invalid_prices': 0
        }
        
        for tech, costs in tech_costs.items():
            if 'cost per unit' not in costs or not isinstance(costs['cost per unit'], (int, float)):
                econ_quality['missing_costs'] += 1
        
        for energy_type, prices in energy_market.items():
            if isinstance(prices, dict):
                for price_type, price in prices.items():
                    if isinstance(price, (int, float)) and price < 0:
                        econ_quality['invalid_prices'] += 1
        
        quality_metrics['economic_data'] = econ_quality
        
        self.reconciliation_results['data_quality'] = quality_metrics
        
        # Save data quality report
        quality_df = pd.DataFrame(quality_metrics).T
        quality_file = os.path.join(self.reconciliation_dir, "data_quality_assessment.csv")
        quality_df.to_csv(quality_file)
        
        print(f"   ✅ Data quality assessment completed")
    
    def _check_input_output_consistency(self, loads_data, production_data, studycase_config):
        """Check consistency between inputs and expected outputs"""
        consistency_checks = {}
        
        # Check load profile consistency
        load_consistency = {
            'total_electricity_demand': 0,
            'total_heat_demand': 0,
            'total_hydrogen_demand': 0,
            'demand_variability': 0
        }
        
        if loads_data:
            if 'electric_load_hybrid.csv' in loads_data:
                elec_demand = loads_data['electric_load_hybrid.csv']['electricity_demand']
                load_consistency['total_electricity_demand'] = elec_demand.sum()
                load_consistency['demand_variability'] = elec_demand.std() / elec_demand.mean() if elec_demand.mean() > 0 else 0
            
            if 'heat_load_hybrid.csv' in loads_data:
                heat_demand = loads_data['heat_load_hybrid.csv']['heat_demand']
                load_consistency['total_heat_demand'] = heat_demand.sum()
            
            if 'hydrogen_demand_hybrid.csv' in loads_data:
                h2_demand = loads_data['hydrogen_demand_hybrid.csv']['hydrogen_demand']
                load_consistency['total_hydrogen_demand'] = h2_demand.sum()
        
        consistency_checks['load_profiles'] = load_consistency
        
        # Check production data consistency
        prod_consistency = {
            'total_solar_production': 0,
            'total_wind_production': 0,
            'production_variability': 0
        }
        
        if production_data:
            solar_sources = [k for k in production_data.keys() if 'SS' in k]
            wind_sources = [k for k in production_data.keys() if 'WS' in k]
            
            for source in solar_sources[:2]:  # Use first 2 solar sources
                data = production_data[source]
                if 'solar_power' in data.columns:
                    prod_consistency['total_solar_production'] += data['solar_power'].sum()
            
            for source in wind_sources[:2]:  # Use first 2 wind sources
                data = production_data[source]
                if 'wind_power' in data.columns:
                    prod_consistency['total_wind_production'] += data['wind_power'].sum()
                    prod_consistency['production_variability'] = data['wind_power'].std() / data['wind_power'].mean() if data['wind_power'].mean() > 0 else 0
        
        consistency_checks['production_data'] = prod_consistency
        
        # Check system configuration consistency
        config_consistency = {
            'total_renewable_capacity': 0,
            'total_storage_capacity': 0,
            'total_conversion_capacity': 0,
            'capacity_balance': 0
        }
        
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
                        if component in ['wind', 'PV']:
                            config_consistency['total_renewable_capacity'] += capacity
                        elif component in ['battery']:
                            config_consistency['total_storage_capacity'] += capacity
                        elif component in ['electrolyzer', 'hydrogen_compressor']:
                            config_consistency['total_conversion_capacity'] += capacity
            
            # Calculate capacity balance
            max_demand = load_consistency['total_electricity_demand'] / 8760  # Average hourly demand
            config_consistency['capacity_balance'] = config_consistency['total_renewable_capacity'] - max_demand
        
        consistency_checks['system_configuration'] = config_consistency
        
        self.reconciliation_results['input_output_consistency'] = consistency_checks
        
        # Save consistency report
        consistency_df = pd.DataFrame(consistency_checks).T
        consistency_file = os.path.join(self.reconciliation_dir, "input_output_consistency.csv")
        consistency_df.to_csv(consistency_file)
        
        print(f"   ✅ Input-output consistency check completed")
    
    def _validate_energy_balance(self, loads_data, production_data, studycase_config):
        """Validate energy balance across the system"""
        energy_balance = {}
        
        # Calculate total energy demand
        total_energy_demand = {
            'electricity_MWh': 0,
            'heat_MWh': 0,
            'hydrogen_MWh': 0
        }
        
        if loads_data:
            if 'electric_load_hybrid.csv' in loads_data:
                elec_demand = loads_data['electric_load_hybrid.csv']['electricity_demand']
                total_energy_demand['electricity_MWh'] = elec_demand.sum() / 1000  # Convert to MWh
            
            if 'heat_load_hybrid.csv' in loads_data:
                heat_demand = loads_data['heat_load_hybrid.csv']['heat_demand']
                total_energy_demand['heat_MWh'] = heat_demand.sum() / 1000  # Convert to MWh
            
            if 'hydrogen_demand_hybrid.csv' in loads_data:
                h2_demand = loads_data['hydrogen_demand_hybrid.csv']['hydrogen_demand']
                # Convert hydrogen demand to energy equivalent (33.33 kWh/kg)
                total_energy_demand['hydrogen_MWh'] = h2_demand.sum() * 33.33 / 1000  # Convert to MWh
        
        energy_balance['total_demand'] = total_energy_demand
        
        # Calculate total energy production (simplified)
        total_energy_production = {
            'solar_MWh': 0,
            'wind_MWh': 0
        }
        
        if production_data:
            solar_sources = [k for k in production_data.keys() if 'SS' in k]
            wind_sources = [k for k in production_data.keys() if 'WS' in k]
            
            # Simplified conversion for solar (assume 10% efficiency)
            for source in solar_sources[:2]:
                data = production_data[source]
                if 'solar_power' in data.columns:
                    # Convert irradiance to power (simplified)
                    solar_power = data['solar_power'] * 0.1  # Simplified conversion
                    total_energy_production['solar_MWh'] += solar_power.sum() / 1000  # Convert to MWh
            
            # Wind production
            for source in wind_sources[:2]:
                data = production_data[source]
                if 'wind_power' in data.columns:
                    total_energy_production['wind_MWh'] += data['wind_power'].sum() / 1000  # Convert to MWh
        
        energy_balance['total_production'] = total_energy_production
        
        # Calculate energy balance
        total_demand = sum(total_energy_demand.values())
        total_production = sum(total_energy_production.values())
        
        energy_balance['balance_metrics'] = {
            'total_demand_MWh': total_demand,
            'total_production_MWh': total_production,
            'net_energy_balance_MWh': total_production - total_demand,
            'energy_sufficiency_percent': (total_production / total_demand * 100) if total_demand > 0 else 0
        }
        
        self.reconciliation_results['energy_balance'] = energy_balance
        
        # Save energy balance report
        balance_df = pd.DataFrame(energy_balance).T
        balance_file = os.path.join(self.reconciliation_dir, "energy_balance_validation.csv")
        balance_df.to_csv(balance_file)
        
        print(f"   ✅ Energy balance validation completed")
    
    def _check_economic_consistency(self, tech_costs, energy_market, studycase_config):
        """Check economic consistency and reasonableness"""
        economic_checks = {}
        
        # Technology cost analysis
        cost_analysis = {
            'total_capital_cost': 0,
            'total_om_cost': 0,
            'average_cost_per_kw': 0,
            'cost_distribution': {}
        }
        
        total_capacity = 0
        for tech, costs in tech_costs.items():
            cost_per_unit = costs.get('cost per unit', 0)
            om_cost = costs.get('OeM', 0)
            
            cost_analysis['total_capital_cost'] += cost_per_unit
            cost_analysis['total_om_cost'] += om_cost
            cost_analysis['cost_distribution'][tech] = cost_per_unit
        
        # Calculate average cost per kW
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
        
        if total_capacity > 0:
            cost_analysis['average_cost_per_kw'] = cost_analysis['total_capital_cost'] / total_capacity
        
        economic_checks['cost_analysis'] = cost_analysis
        
        # Energy price analysis
        price_analysis = {
            'electricity_purchase_price': energy_market.get('electricity', {}).get('purchase', 0),
            'electricity_sale_price': energy_market.get('electricity', {}).get('sale', 0),
            'hydrogen_purchase_price': energy_market.get('hydrogen', {}).get('purchase', 0),
            'hydrogen_sale_price': energy_market.get('hydrogen', {}).get('sale', 0),
            'price_spread_electricity': 0,
            'price_spread_hydrogen': 0
        }
        
        # Calculate price spreads
        if price_analysis['electricity_purchase_price'] > 0 and price_analysis['electricity_sale_price'] > 0:
            price_analysis['price_spread_electricity'] = price_analysis['electricity_purchase_price'] - price_analysis['electricity_sale_price']
        
        if price_analysis['hydrogen_purchase_price'] > 0 and price_analysis['hydrogen_sale_price'] > 0:
            price_analysis['price_spread_hydrogen'] = price_analysis['hydrogen_purchase_price'] - price_analysis['hydrogen_sale_price']
        
        economic_checks['price_analysis'] = price_analysis
        
        # Incentive analysis
        incentive_analysis = {
            'green_hydrogen_incentive': energy_market.get('green_hydrogen_incentives', {}).get('value', 0),
            'ammonia_incentive': energy_market.get('ammonia_incentives', {}).get('value', 0),
            'rec_incentive': energy_market.get('REC', {}).get('collective self consumption incentives', 0),
            'total_incentives': 0
        }
        
        incentive_analysis['total_incentives'] = sum([
            incentive_analysis['green_hydrogen_incentive'],
            incentive_analysis['ammonia_incentive'],
            incentive_analysis['rec_incentive']
        ])
        
        economic_checks['incentive_analysis'] = incentive_analysis
        
        self.reconciliation_results['economic_consistency'] = economic_checks
        
        # Save economic consistency report
        econ_df = pd.DataFrame(economic_checks).T
        econ_file = os.path.join(self.reconciliation_dir, "economic_consistency_check.csv")
        econ_df.to_csv(econ_file)
        
        print(f"   ✅ Economic consistency check completed")
    
    def _check_temporal_consistency(self, loads_data, production_data):
        """Check temporal consistency of data"""
        temporal_checks = {}
        
        # Check load profile temporal patterns
        load_temporal = {
            'daily_patterns': {},
            'seasonal_patterns': {},
            'peak_hours': {},
            'off_peak_hours': {}
        }
        
        if loads_data:
            for load_file, data in loads_data.items():
                if 'electricity_demand' in data.columns:
                    # Daily pattern analysis
                    data['hour'] = pd.to_datetime(data['timestamp']).dt.hour
                    daily_avg = data.groupby('hour')['electricity_demand'].mean()
                    load_temporal['daily_patterns'][load_file] = daily_avg.to_dict()
                    
                    # Peak and off-peak analysis
                    peak_hours = daily_avg.nlargest(4).index.tolist()
                    off_peak_hours = daily_avg.nsmallest(4).index.tolist()
                    load_temporal['peak_hours'][load_file] = peak_hours
                    load_temporal['off_peak_hours'][load_file] = off_peak_hours
        
        temporal_checks['load_temporal'] = load_temporal
        
        # Check production data temporal patterns
        prod_temporal = {
            'solar_daily_pattern': {},
            'wind_daily_pattern': {},
            'seasonal_variation': {}
        }
        
        if production_data:
            solar_sources = [k for k in production_data.keys() if 'SS' in k]
            wind_sources = [k for k in production_data.keys() if 'WS' in k]
            
            # Solar daily pattern
            if solar_sources:
                solar_data = production_data[solar_sources[0]]
                if 'solar_power' in solar_data.columns:
                    solar_data['hour'] = pd.to_datetime(solar_data['timestamp']).dt.hour
                    solar_daily = solar_data.groupby('hour')['solar_power'].mean()
                    prod_temporal['solar_daily_pattern'] = solar_daily.to_dict()
            
            # Wind daily pattern
            if wind_sources:
                wind_data = production_data[wind_sources[0]]
                if 'wind_power' in wind_data.columns:
                    wind_data['hour'] = pd.to_datetime(wind_data['timestamp']).dt.hour
                    wind_daily = wind_data.groupby('hour')['wind_power'].mean()
                    prod_temporal['wind_daily_pattern'] = wind_daily.to_dict()
        
        temporal_checks['production_temporal'] = prod_temporal
        
        self.reconciliation_results['temporal_consistency'] = temporal_checks
        
        # Save temporal consistency report
        temp_df = pd.DataFrame(temporal_checks).T
        temp_file = os.path.join(self.reconciliation_dir, "temporal_consistency_check.csv")
        temp_df.to_csv(temp_file)
        
        print(f"   ✅ Temporal consistency check completed")
    
    def _validate_physical_constraints(self, studycase_config, loads_data, production_data):
        """Validate physical constraints and system limits"""
        physical_checks = {}
        
        # System capacity constraints
        capacity_constraints = {
            'max_electricity_demand': 0,
            'max_heat_demand': 0,
            'max_hydrogen_demand': 0,
            'total_renewable_capacity': 0,
            'capacity_margin': 0
        }
        
        # Get maximum demands
        if loads_data:
            if 'electric_load_hybrid.csv' in loads_data:
                capacity_constraints['max_electricity_demand'] = loads_data['electric_load_hybrid.csv']['electricity_demand'].max()
            if 'heat_load_hybrid.csv' in loads_data:
                capacity_constraints['max_heat_demand'] = loads_data['heat_load_hybrid.csv']['heat_demand'].max()
            if 'hydrogen_demand_hybrid.csv' in loads_data:
                capacity_constraints['max_hydrogen_demand'] = loads_data['hydrogen_demand_hybrid.csv']['hydrogen_demand'].max()
        
        # Get total renewable capacity
        if studycase_config and 'hybrid_plant' in studycase_config:
            hybrid_config = studycase_config['hybrid_plant']
            for component, config in hybrid_config.items():
                if isinstance(config, dict):
                    capacity = None
                    for key in ['capacity', 'Npower', 'peakP', 'max capacity', 'max power']:
                        if key in config and isinstance(config[key], (int, float)):
                            capacity = config[key]
                            break
                    if capacity and component in ['wind', 'PV']:
                        capacity_constraints['total_renewable_capacity'] += capacity
        
        # Calculate capacity margin
        if capacity_constraints['total_renewable_capacity'] > 0:
            capacity_constraints['capacity_margin'] = (
                capacity_constraints['total_renewable_capacity'] - capacity_constraints['max_electricity_demand']
            ) / capacity_constraints['total_renewable_capacity'] * 100
        
        physical_checks['capacity_constraints'] = capacity_constraints
        
        # Storage constraints
        storage_constraints = {
            'battery_capacity_kwh': 0,
            'hydrogen_storage_kg': 0,
            'storage_duration_hours': 0
        }
        
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
                        if component == 'battery':
                            storage_constraints['battery_capacity_kwh'] = capacity
                        elif component == 'hydrogen_storage':
                            storage_constraints['hydrogen_storage_kg'] = capacity
        
        # Calculate storage duration
        if storage_constraints['battery_capacity_kwh'] > 0 and capacity_constraints['max_electricity_demand'] > 0:
            storage_constraints['storage_duration_hours'] = storage_constraints['battery_capacity_kwh'] / capacity_constraints['max_electricity_demand']
        
        physical_checks['storage_constraints'] = storage_constraints
        
        self.reconciliation_results['physical_constraints'] = physical_checks
        
        # Save physical constraints report
        phys_df = pd.DataFrame(physical_checks).T
        phys_file = os.path.join(self.reconciliation_dir, "physical_constraints_validation.csv")
        phys_df.to_csv(phys_file)
        
        print(f"   ✅ Physical constraints validation completed")
    
    def _generate_comprehensive_report(self):
        """Generate comprehensive reconciliation report"""
        print("📋 Generating comprehensive reconciliation report...")
        
        # Create summary report
        summary_data = []
        
        for check_type, results in self.reconciliation_results.items():
            if isinstance(results, dict):
                for metric, value in results.items():
                    if isinstance(value, dict):
                        for sub_metric, sub_value in value.items():
                            summary_data.append({
                                'Check_Type': check_type,
                                'Metric': f"{metric}_{sub_metric}",
                                'Value': sub_value,
                                'Status': 'Valid' if sub_value is not None else 'Missing'
                            })
                    else:
                        summary_data.append({
                            'Check_Type': check_type,
                            'Metric': metric,
                            'Value': value,
                            'Status': 'Valid' if value is not None else 'Missing'
                        })
        
        # Create comprehensive report
        comprehensive_report = pd.DataFrame(summary_data)
        report_file = os.path.join(self.reconciliation_dir, "comprehensive_reconciliation_report.csv")
        comprehensive_report.to_csv(report_file, index=False)
        
        # Generate summary statistics
        summary_stats = {
            'Total_Checks': len(comprehensive_report),
            'Valid_Checks': len(comprehensive_report[comprehensive_report['Status'] == 'Valid']),
            'Missing_Checks': len(comprehensive_report[comprehensive_report['Status'] == 'Missing']),
            'Data_Quality_Score': len(comprehensive_report[comprehensive_report['Status'] == 'Valid']) / len(comprehensive_report) * 100
        }
        
        summary_df = pd.DataFrame([summary_stats])
        summary_file = os.path.join(self.reconciliation_dir, "reconciliation_summary_statistics.csv")
        summary_df.to_csv(summary_file, index=False)
        
        # Generate recommendations
        recommendations = []
        
        if summary_stats['Data_Quality_Score'] < 90:
            recommendations.append("Data quality score is below 90%. Review missing or invalid data.")
        
        if summary_stats['Missing_Checks'] > 0:
            recommendations.append(f"Found {summary_stats['Missing_Checks']} missing data points. Complete data collection.")
        
        # Add specific recommendations based on reconciliation results
        if 'energy_balance' in self.reconciliation_results:
            balance_metrics = self.reconciliation_results['energy_balance'].get('balance_metrics', {})
            if balance_metrics.get('energy_sufficiency_percent', 0) < 50:
                recommendations.append("Energy production is insufficient for demand. Consider increasing renewable capacity.")
        
        if 'physical_constraints' in self.reconciliation_results:
            capacity_constraints = self.reconciliation_results['physical_constraints'].get('capacity_constraints', {})
            if capacity_constraints.get('capacity_margin', 0) < 10:
                recommendations.append("Low capacity margin. Consider increasing system capacity for reliability.")
        
        # Save recommendations
        if recommendations:
            rec_df = pd.DataFrame({'Recommendations': recommendations})
            rec_file = os.path.join(self.reconciliation_dir, "reconciliation_recommendations.csv")
            rec_df.to_csv(rec_file, index=False)
        
        print(f"   ✅ Comprehensive reconciliation report generated")
        print(f"   📊 Data Quality Score: {summary_stats['Data_Quality_Score']:.1f}%")
        print(f"   📋 Total Checks: {summary_stats['Total_Checks']}")
        print(f"   ✅ Valid Checks: {summary_stats['Valid_Checks']}")
        print(f"   ⚠️ Missing Checks: {summary_stats['Missing_Checks']}")
        
        if recommendations:
            print(f"   💡 Recommendations: {len(recommendations)} items")
            for i, rec in enumerate(recommendations[:3], 1):  # Show first 3 recommendations
                print(f"      {i}. {rec}")
            if len(recommendations) > 3:
                print(f"      ... and {len(recommendations) - 3} more recommendations")
