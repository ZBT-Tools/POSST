# -*- coding: utf-8 -*-
"""
**Python Open Source Simulation Tool (POSST): Simulating Transient Systems Based on Controlled-Oriented Models (COM), Mean Value Modeling (MVM), Lumped Parameter System (LPS)**

### Authors:
Florian Dennewitz (Maintainer), Leonora Kastrati (Maintainer), Niklas Nickig (Maintainer), Matthias Bahr (Contributor)

---

## Abstract
This open-source simulation software enables real-time simulation of fuel cell system components without the use of mathematical solvers. Based on the Mean-Value-Modelling approach for black box modelling, this software is web-based and addresses common challenges of existing simulation programs, such as high costs and lack of interoperability.

### Example of Cathode Path

<p align="center">
  <img src="https://zbt-tools.github.io/POSST/pic/figure1.png" alt="Overview of Cathode Path" width="500">
</p>

### Results of Example
<p align="center">
  <img src="https://zbt-tools.github.io/POSST/pic/figure8.png" alt="Cathode Load Profile Simulation" width="800">
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
  <img src="https://zbt-tools.github.io/POSST/pic/figure2.png" alt="Reservoir Inputs and Outputs" width="500">
</p>

### 2. Modelling of Flows
Assuming incompressible flow between reservoirs, mass flow is calculated with adaptations for potential backflow.

### 3. Modelling of Mixed Media
To calculate specific heat capacity and other properties, the gas composition is determined.

### 4. Modelling of the Compressor
The compressor is simulated using the “rotating disk” model. A Python library generates the compressor characteristic diagram.

<p align="center">
  <img src="https://zbt-tools.github.io/POSST/pic/figure3.png" alt="Compressor Model - Rotating Disk Approach" width="500">
</p>

A Python library generates the compressor characteristic diagram:

<p align="center">
  <img src="https://zbt-tools.github.io/POSST/pic/figure5.png" alt="Compressor Characteristic Diagram" width="500">
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
- [Hydrogen and Fuel Cell Center (ZBT), Duisburg, Germany](mailto:f.dennewitz@zbt.de)
- [Energy System Solutions (ESS), Bottrop, Germany](mailto:florian.dennewitz@unitybox.de)

---

"""
