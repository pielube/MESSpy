# MESSpy - Multi-Energy System Simulator - Python version
### Authors
**Pietro Lubello** Università degli Studi di Firenze, Italy\
**Mattia Pasqui** Università degli Studi di Firenze, Italy\
**Alessandro Mati** Università degli Studi di Firenze, Italy\
**Carlo Carcasci** Università degli Studi di Firenze, Italy

### Overview
Multi-Energy System Simulator born to make techno-economic assesment of Renewable Energy Communities (REC), but can also be used to study a single building.
It simulate hour after hour the energy flows between technologies of each location (building) inside the REC and calculate the Net Present Value of each considering relationship with the national grid and incentives. The program has been developed to be as generalizable as possible so it can be used to simulate even very different case studies and easily change their configuration or characteristics, technical and economic.
The program is fully commented and can be easly used as black box by modifing input/ and working on results/ or also modifyng the code.

It's an objected oriented program structured on three levels: REC, location and technologies.
Models of different technologies are aviable and still under developing to include new feture and more realistic detalis. At the moment the following technologies can be included in the simulations:
- Photovoltaic field
- Battery
- Electrolyzer
- Fuel Cells
- Hydrogen Tank
- Heat Pumps

### Requirements
The model is developed in Python 3.9, and requires the following libraries:
- numpy
- pandas
- os
- pickle (results are saved in .pickle)
- json (input files are .json)

The following libraries are useful but not necessary, so you can not use them by giving up some functionality:
- time (used only to check code speed)
- pvlib (used to download PV production series based on typical meteorological year)
- matplotlib (used in post_process)
- plotly (used in post_process.flow())

### Quick start
To get started, download the repository and simply run the "run_test.py" script

### Input files
You can modofy them from a python interface or simply from notepad. Input_test/ contains contains a demonstration case study. 
- general.json defines the general input. More details can be found in rec.py following the early comments
- structure.json defines the structure of the case study. Here you can define all the locations to consider, each technologie inside the locations and technology's parameters. More detalis can be found following comments in rec.py and location.py.
- refcase.json This file has the same structure of structure.json and defines the "buiseness as usual" case, which is used as a reference case to calculate cash flows of the study case and make economis assesment.
- economics.json defines economic parameters. More details can be found following comments in economics.py

### How to continue
We suggest to create your own run_dev.py, input_dev/ and post_process_dev.py and work on them instead of modify the existing file using as initial test. 

### Citing
Please cite the original Journal publication if you use MESSpy in your research: 
work in progress ...

### Related works
- Renewable Energy Communities: a techno-economic assesment focusing on smart battery management
- Work in progres ..

### Contribute
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
Copyright 2022 MESSpy, contributors listed in Authors

Licensed under the European Union Public Licence (EUPL), Version 1.2-or-later; you may not use this file except in compliance with the License.

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License
