# POSST
**Python Open Source Simulation Tool (POSST): Simulating Transient Systems Based on Controlled-Oriented Models (COM), Mean Value Modeling (MVM), Lumped Parameter System (LPS)**

### Authors:
Florian Dennewitz (Maintainer), Leonora Kastrati (Maintainer), Niklas Nickig (Maintainer), Matthias Bahr (Contributor)

---

## Abstract
This open-source simulation software enables real-time simulation of fuel cell system components without the use of mathematical solvers. Based on the Mean-Value-Modelling approach for black box modelling, this software is web-based and addresses common challenges of existing simulation programs, such as high costs and lack of interoperability.

### Example of Cathode Path
<p align="center">
  <img src="pic/figure1.png" alt="Overview of Cathode Path" width="500">
</p>

### Results of Example
<p align="center">
  <img src="pic/figure8.png" alt="Cathode Load Profile Simulation" width="800">
</p>

## Keywords
- Simulation
- Open-source
- Python
- Black box modelling
- Fuel cell system

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
1. Zahn, S. “Arbeitsspielaufgelöste Modellbildung und Hardware-in-the-Loop-Simulation von Pkw-Dieselmotoren mit Abgasturboaufladung.”
2. Guzzella, L., Onder, C. “Introduction to Modeling and Control of Internal Combustion Engine Systems.”
3. Virtanen, P. et al. “SciPy 1.0: Fundamental Algorithms for Scientific Computing in Python.”
4. Gößling, S. “2-D + 1-D ortsaufgelöste Modellierung von PEM-Brennstoffzellen.”

---

### Contact
- [Hydrogen and Fuel Cell Center (ZBT), Duisburg, Germany](mailto:f.dennewitz@zbt.de)
- [Energy System Solutions (ESS), Bottrop, Germany](mailto:florian.dennewitz@unitybox.de)
