# MESSpy - Multi-Energy System Simulator
An open-source model for simulating multy-energy systems such as Renewable Energy Communities and much more.

![image](https://user-images.githubusercontent.com/83342584/177154238-c0e6b9fa-69b4-40c8-9b45-9d16919e079b.png)

### Authors
**Pietro Lubello**, **Mattia Pasqui**, **Alessandro Mati**, **Andrea Ademollo** and **Carlo Carcasci**\
Department of Industrial Engineering (DIEF), University of Florence (UNIFI), Italy

## Overview
The Multi-Energy System Simulator has been developed to perform techno-economic assesment of Renewable Energy Communities (REC), but can also be used to study single-standing buildings and hydrogen-integrated energy systems.
It can simulate hourly balances of the energy flows between technologies of each location (building) inside the REC and calculate the Net Present Value of each considering the interaction with the national grid given different incentive schemes. The program has been developed to be as general as possible so it can be used to simulate a wide range of different case studies while easily changing their configuration or parameters,both technical and economic.
The code is extensively commented and can be easily used either as a black box by simply modifying the inputs and working on results or by directly modifying the code.

It's an object oriented program structured on three levels: REC, location and technologies.
Models of different technologies are avialable and still under development to include new fetures and more realistic detalis. At the current stage the following technologies can be included in the simulations:
- Photovoltaic panels
- Wind turbines
- Batteries
- Electrolyzers
- Fuel Cells
- Hydrogen tanks
- Heat pumps
- Boilers (ng, ele or H2)
- CHP (combained heat and power)

### MESS needs the hourly load profiles as input as a .csv file
Depending on the type of meter installed, these data is in some cases made available by the electricity/gas supplier, in others it must be requested, while sometimes it cannot be obtained. In the latter case, specific programmes are required to generate such profiles in the specific .csv format needed as one of the program inputs. There are many programmes available online, the authors recommend the following:\
bottom-up model: https://github.com/RAMP-project/RAMP \
top-down model: https://github.com/PasquinoFI/LoBi

### Requirements
The model is developed in Python 3.9, and requires the following libraries:
- numpy
- pandas
- os
- pickle (results are saved in .pickle)
- json (input files are .json)
- pvlib (used to download PV production series and weather data based on typical meteorological year)
- matplotlib (used in post_process)

A less up-to-date but fully functional and documented fortran version is also available:
https://github.com/pielube/MESS-Fortran

## Quick start
To get started, download the repository and simply run the "run_test.py" script

### Input files
You can modify them from a python interface or simply from notepad. The "Input_test" folder contains a demonstration case study. 
- general.json defines the general input. More details can be found in rec.py prologue comments.
- structure.json defines the structure of the case study. Here you can define all the locations to consider, each technology inside the locations and technology's parameters. More detalis can be found in rec.py and location.py comments to the code.
- refcase.json This file has the same structure of structure.json and defines the "buiseness as usual" case, which is used as a reference case for calculating the cash flows of the study case and performing the economic assessment.
- economics.json defines economic parameters. More details can be found in the comments of economics.py

### How to continue
We suggest you to create your own run_dev.py, input_dev/ and post_process_dev.py and to work on them instead of modifying the existing file used as initial test. 

## Related works
- "Optimal sizing of a distributed energy system with thermal load electrification"\
https://www.e3s-conferences.org/articles/e3sconf/abs/2020/57/e3sconf_ati2020_01006/e3sconf_ati2020_01006.html
- "The potential of simulating energy systems: The multi energy systems simulator model"\
https://www.mdpi.com/1996-1073/14/18/5724
- "Considerations on the impact of battery ageing estimation in the optimal sizing of solar home battery systems" https://www.sciencedirect.com/science/article/pii/S0959652621039299
- "Assessment of hydrogen-based long term electrical energy storage in residential energy systems"
https://www.sciencedirect.com/science/article/pii/S2666955222000260
- "A New Smart Batteries Management for Renewable Energy Communities"                                                   
https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4268979
- "Decentralised heat pumps management in Renewable Energy Communities"

## Citing
Please cite previous works if you use MESSpy in your research.

## Contribute
This project is open-source. Interested users are therefore invited to test, comment or contribute to the tool. Submitting issues is the best way to get in touch with the development team, which will address your comment, question, or development request in the best possible way. We are also looking for contributors to the main code, willing to contibute to its capabilities, computational-efficiency, formulation, etc.

To contribute changes:

- Fork the project on GitHub
- Create a feature branch (e.g. named "add-this-new-feature") to work on in your fork
- Commit your changes to the feature branch
- Push the branch to GitHub
- On GitHub, create a new pull request from the feature branch
- When committing new changes, please also take care of checking code stability running run_test 
- Your name will be added to the authors list

### License
Copyright 2022 MESSpy, contributors listed in Authors.

Licensed under the European Union Public Licence (EUPL), Version 1.2-or-later; you may not use this file except in compliance with the License.

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.
