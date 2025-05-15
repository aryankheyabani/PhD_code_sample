# PhD Code Sample 
This is a sample of the in-house MATLAB code which I developed during my PhD.

The purpose of this repository is to show my skills in coding in MATLAB for numerical modeling and computational mechanics. 

The code included here is a part of the script I developed to implement the Parametric High Fidelity Generalized Method of Cells (PHFGMC) micromechanics. 

The PHFGMC is a micromechanics method to obtain the macroscale homogenized properties of a periodic composite material based on the properties of its constituent materials at the microscale. To compute the macroscale material properties, a repeating unit cell (RUC) is generated based on the periodic structure of the composite material. Then, the required computations are performed on the RUC to get the material properties. The RUC demonstrates how the constituent materials (fiber and matrix in this case) form the microstructure, geometries of those material as well as the fiber volume fraction (amount of fiber per total volume of the RUC). So for the PHFGMC, a geometry model is needed to define the microstructure. This model is then discretized into subcells similar to elements of a FEM RVE (representative volume element) to include the different constituent materials. 

The PHGMC formulation applies different boundary conditions on the RUC to get the system of governing equations. These include continuity of displacements and tractions between the subcells and equilibrium equations in each subcell. Displacement and tractions continuities are applied in an average integral sense between the faces of subcells inside the RUC. To this end, the code needs those face connectivities as input. Herein, the geometry of the RUC is modeled and the face connectivities are already defined and given as input to the code. Furthermore, material properties for the constituent materials need to be provided.

Note: This code is not functional and not intended to work since all the necessary functions are not included. 
