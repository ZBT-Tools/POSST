


# Python Open Source Simulation Tool (POSST):

**Simulating Transient Systems Based on Controlled-Oriented Models (COM), Mean Value Modeling (MVM), Lumped Parameter System (LPS)**

### Project partner:

<p align="left">
  <img src="pic/ESS_LOGO_vektor.svg" alt="ESS Logo" width="40">
<a href="https://www.energy-system-solutions.de/">Energy System Solutions (ESS)<sup>1</sup></a>
</p>
<p align="left">
  <img src="pic/ZBT_Logo_RGB_B_L.svg" alt="ZBT Logo" width="40">
<a href="https://zbt.de/das-zbt/wissenschaftliche-abteilungen/brennstoffzellensysteme/">Hydrogen and Fuel Cell Center (ZBT)<sup>2</sup></a>
</p>

[Become a project partner](mailto:florian.dennewitz@unitybox.de)

### Authors:
Florian Dennewitz (Maintainer)<sup>1,2</sup>, Leonora Kastrati (Maintainer)<sup>1</sup>, Niklas Nickig (Maintainer)<sup>2</sup>, Matthias Bahr (Contributor)<sup>2</sup>

---

## Abstract
This open-source simulation software enables real-time simulation of fuel cell system components without the use of mathematical solvers. Based on the Mean-Value-Modelling approach for transient process modelling, this software is web-based and addresses common challenges of existing simulation programs, such as high costs and lack of interoperability.

### Example of Cathode Path
<p align="center">
  <img src="pic/figure1.png" alt="Overview of Cathode Path" width="500">
</p>

### Results of Example
<p align="center">
  <img src="pic/figure8.png" alt="Cathode Load Profile Simulation" width="800">
</p>

### Link to Websimulation
Here we provide a simulation environment on the Energy System Solutions Website:
https://www.energy-system-solutions.de/simulation. Feel free to test the simulation and the integrated controller of the ZBT.

---

## Keywords

- Simulation
- Open-source
- Python
- Black box modelling
- Fuel cell system

---

## Installation Guide
### Installation with pip
NOT AVAILABLE YET!


### Installation with GIT

If you prefer to install directly from the repository, follow these steps:

Clone the Repository:

    git clone https://github.com/ZBT-Tools/POSST.git

Navigate to the Project Directory:

    cd yourPath/POSST

Install Required Packages: Make sure you have a virtual environment set up (optional but recommended). Then install the required packages:

    pip install -r requirements.txt

Create Your Process: Follow the documentation in the repository to set up your process model.

Simulate Process: Use the provided simulation tools to run your process model. Check the usage examples in the documentation for guidance.

Additional Notes

    Ensure that you have Python installed (preferably Python 3.7 or higher).
    If you encounter any issues, check the issues section of the repository for troubleshooting.

Feel free to reach out if you have any questions or need further assistance!

---

## Contents

### 1. Modelling of Reservoirs
The reservoir model describes mass and enthalpy flow as well as temperature changes using ordinary differential equations.

<p align="center">
  <img src="pic/figure2.png" alt="Reservoir Inputs and Outputs" width="500">
</p>

### 2. Modelling of Flows
Assuming incompressible flow between reservoirs, mass flow is calculated with adaptations for potential backflow.

### 3. Modelling of Mixed Media
To calculate specific heat capacity and other properties, the gas composition is determined.

### 4. Modelling of the Compressor
The compressor is simulated using the “rotating disk” model. A Python library generates the compressor characteristic diagram.

<p align="center">
  <img src="pic/figure3.png" alt="Compressor Model - Rotating Disk Approach" width="500">
</p>

A Python library generates the compressor characteristic diagram:

<p align="center">
  <img src="pic/figure5.png" alt="Compressor Characteristic Diagram" width="500">
</p>

### 5. Modelling of the Fuel Cell
To reduce a 2-D + 1-D spatially resolved fuel cell model to a voltage model, several simplifications are necessary. By focusing only on the voltage model and ignoring certain sub-models, the complexity of the simulation can be significantly reduced while maintaining core performance characteristics relevant to voltage output.

---

## References
[1]	Dipl.-Ing. Sebastian Zahn, "Arbeitsspielaufgelöste Modellbildung und Hardware-in-the-Loop-Simulation von Pkw-Dieselmotoren mit Abgasturboaufladung," Fachbereich Elektrotechnik und Informationstechnik, Technische Universität Darmstadt, Darmstadt, 2012.<br>
[2]	L. Guzzella und C. Onder, Introduction to Modeling and Control of Internal Combustion Engine Systems, 2. Aufl. Berlin, Heidelberg: Springer Berlin Heidelberg, 2009. Zugriff am: 29. Oktober 2024. [Online]. Verfügbar unter: https://link.springer.com/book/10.1007/978-3-642-10775-7<br>
[3]	P. Virtanen et al., "SciPy 1.0: Fundamental Algorithms for Scientific Computing in Python," Nature Methods, Jg. 17, S. 261–272, 2020, doi: 10.1038/s41592-019-0686-2.<br>
[4]	S. Gößling, "2-D + 1-D ortsaufgelöste Modellierung von PEM-Brennstoffzellen," Fakultät für Ingeneuirwissenschaften, Abteilung Maschienenbau und Verfahrenstechnik, Universität Duisburg Essen, Duisburg, 2019.<br>



---

### Contact
- [Energy System Solutions (ESS), Bottrop, Germany](mailto:florian.dennewitz@unitybox.de)
- [Hydrogen and Fuel Cell Center (ZBT), Duisburg, Germany](mailto:f.dennewitz@zbt.de)

## Documentation
- [Docu](https://zbt-tools.github.io/POSST/html/POSST.html)