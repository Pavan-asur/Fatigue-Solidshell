# Fatigue-Solidshell
The codes presented here are used for the following paper "A Phase field fracture model for fatigue using Locking free solid shell finite elemnts: Analysis for homogenous and layered composites".
The codes can be used with or without fatigue. The code is written as a UEL and UMAT for the use in abaqus. We have tested the files for version Abaqus2020, and Abaqus2022. 

Do not forget to change the "allelem" with number of the elements in the UEL in each layer. The variable   "allelem" is located in the module Kvisual. 

Care has to be taken when defining the fatigue model (from the input file).

If you are using this code for academic and research purpose. Also check out other codes on Thermo-mechanical solid shell, and do not forget to cite our article. 

Pavan Kumar Asur Vijaya Kumar, Aamir Dean, Jose Reinoso, Heinz Pettermann, Marco Paggi,  
“A phase-field fracture model for fatigue using locking-free solid shell finite elements: Analysis for homogeneous materials and layered composites”. 
In: Theoretical and Applied Fracture Mechanics (2023), p. 104029. issn: 0167-8442. doi: https://doi.org/10.1016/j.tafmec.2023.104029.


Please report any problems and suggestions to pavan.kumar@ilsb.tuwien.ac.at

Co-Authors:

Jose Reinoso, "josereinoso", University of Seville, Spain.

Aamir Dean, "aamir-dean", Leibniz University of Hannover, Germany. 
