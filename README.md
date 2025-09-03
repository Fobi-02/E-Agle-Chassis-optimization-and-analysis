# E-Agle-Chassis-optimization-and-analysis
One of the most important charatteristic of a chassis is it’s torsional stiffness because we want it to be as rigid as possible to make the suspension properly work as intended.

For the Kraken chassis we wanted to develope a tool that could help us to find the best position of the nodes to have the best stiffness possible without having to rely on manually setup many simulation on ansys mechanical.

## For what is this script useful?
This script was used not only to automatically optimize the chassis but also:

* To quickly iterate different chassis design without setting up many ansys simulations
* To analize the forces in all tubes
* To calulate the forces acting on the rocker for NTop topological optimization
* To calculate forces acting on joints in case of an impact to later analize with FEM
* To analize the forces in tubes for press fit sizing

## Working principle
This Wolfram Mathematica script has the goal of calculating numerically the torsional stiffness of a certain chassis configuration by only giving the nodes coordinates and tubes propertuies (start and end position, material, diameter…) as an input.

Using Mathematica optimization method (Nelder Mead in this case) it is possible to use this alghorithm to optimize the location of the nodes to obtain the best stiffness.

This script builds a mathematical model of the structure and applyies some forces to it. Then it solves the model to measure the displacement and so the stiffness. The main steps of this script are:

\-Solving the externally isostatic reaction forces

\-Use node to node analysis to write equation to calculate the interal reaction forces

\-Calculate the remaining reaction forces using the Castigliano method.

\-Use Mohr Maxwell to compute the displacement of the node of interest

\-Calculate the torsional stiffness given the displacement.

### Initial assumptions
We considered the chassis members as Saint Venant beams that are subject to forces only at the start and at the end of them. Since in a formula student chassis is not possible to have the beams work only in tension and compression (not possible to fully triangulate the structure) we have modelled the nodes between beams as rigid joints that fully blocks all the degrees of freedom. For this reason the structure will become many times hyperstatic.

## Explaination example
To better understand the working principle of this model it is useful to look at a simpler 2D example of an internal hyperstatic structure subject to  one force. This triangle made of three deformable beams is fully blocked in A.

<img width="1412" height="410" alt="image" src="https://github.com/user-attachments/assets/4d9f0e20-7fce-4d65-bdd6-980a40e5dbb1" />


Based on the geometry we can directly calculate the reaction forces (since it is externally isostatic).

We can then write equilibrium equation for each beam and each node in the structure.

<img width="1451" height="795" alt="image (1)" src="https://github.com/user-attachments/assets/ca654f0f-b825-4d06-81ce-f704ba197422" />

 As can be seen the N and Ty reaction forces are symmetric along the beam while Mx is different at the start and at the end of a beam. For each beam we can only write the rotational equilibrium (the other equilibrium are trivial) while for each node we can write the two translational equilibrium and the rotational equilibrium.

For this example we have 12 unknown variables, 9 equations for the nodes and 3 equations for the beams. But we can not directly calculate the internal forces since the system is hyperstatic. Trying to solve this system with wolfram mathematica we can find solution only to 9 out of the 12 variables. The other 3 must be solved using castigliano.

I can choose the variables to leave as parameter, in this case i have considered N3, Ty2, Mx1A.


I have to calculate the internal action to be able to apply the Castigliano theorem later on.

<img width="691" height="622" alt="image (2)" src="https://github.com/user-attachments/assets/fe6af6a7-bfeb-49a5-9d1c-882c7e753896" />


Writing a system with 3 castigliano equations i can solve the remaining internal reaction forces:

<img width="527" height="856" alt="image (3)" src="https://github.com/user-attachments/assets/58618e76-fbf1-4d20-a73c-7d9d90707561" />


We can then apply different methods to calculate the displacement of a certain point. For example we can use the castigliano theorem leving the force as a variable and differentiating in respect to it. For this algorithm tho it has been used the Mohr Maxwell theorem to be able to solve the system of equations numerically and not simbolically. For this method we need to apply an unit force where we want to mesure the displacement and solve the internal force of that sistem

## Mathematical solutions

Before explaining the code is important to understand some solutions used in this code

### Component directions
<img width="1242" height="1076" alt="image" src="https://github.com/user-attachments/assets/e96c10f0-7077-4351-a2ac-20d69aa865ef" />

I used the right hand rule for the forces at the end of the beam, the forces at the start of the beams do not follow this rule because they need to be opposite to the one on the end point. For this model we have no forces applied in the middle of a beam (if we need to we just use two beam one after the other) so the module of N, Tx, Ty, Mz are equal from the start point to the end point. The only difference is the Mx and My that changes along the beam since we have the shear forces applied. The moment at the start of the beam can be easily calculated computing the rotational equilibrium. So knowing this for each beam we have 6 unknown variables to calculate.

We need a standard definition for the forces direction. I have decided to locate the z axis along the beam length. Then to find an easily recognizable standart for the y axis i have decided to have:

* y perpendicular to z
* y parallel to the z=0 plane
* Ty component (of the end node)  with y component positive

  The x axis instead is defined as:
* x perpendicular to z
* x perpendicular to y
* x,y,z follows the right hand rule

The only exception to this notation is when the beam is straigt perpendicular to the z=0 plane. In this case we can not have the y axis perpendicular to that plane so we chose as direction the vector {0, 1, 0}

### Computing the direction

Knowing the coordinates of the start node (A) and the end node (B) for a beam i can compute the direction of the beam using a function that i will call dir\[A,B\].

We want to calculate the components of the vector for each force/moment starting from his module.

I will talk about how to calculate the values for the end node of the beam, if we want to calculate the values in the starting node we simply have to multiply by -1 to invert the sign.

<img width="1128" height="1110" alt="image" src="https://github.com/user-attachments/assets/2759b408-7480-4d5d-b82a-da82badb426d" />

We want to calculate all the directions as explained above. We need to calculate the vector V and W with unit mdule. This versor will multiply the module of the forces to calculate their components in the global (x,y,z) coordintes.

The minus sign before all the values is due to the fact that the forces exiting the beam are entering the node.

### Node equations

In order to write the equilibrium equation i need for each node to:

* Identify the connected beams (the information is stored in a table - connection matrix)
* Check if the node is the starting point or the ending point (information stored in the beam matrix)
* sum the forces: N\[k\]+Tx\[k\]+Ty\[k\] with k every beam connected
* Sum the moments: Mx\[k\]+Mx\[k\]+My\[k\]
* Add the external forces acting on the node
* We will have 6 equations for each node since each element is a vector with 3 components 

  (ex. N\[2\]={1, 0, 0})
* All the equations for each node are collected into a system

### Internal actions

Internal actions are quite easy to compute once the vectors of the forces have been defined

<img width="1337" height="776" alt="image" src="https://github.com/user-attachments/assets/01c5fa11-fc46-40af-a8bf-7d41f0843adc" />

### Castigliano

To make the calculation more efficients and to save compute time i have manually solved the integral knowing how Mx and My are defined.

<img width="1382" height="969" alt="image" src="https://github.com/user-attachments/assets/ae567115-06cf-44d6-a8d7-3f7823392f57" />

## Code explaination step by step

### Problem setup
In this section of the code we need to give the system all the informations it needs to build the chassis mathematical model.

### Nodes matrix

This table contains the coordinates and the forces applied to each node of the chassis. The table is composed as follows:

(nodeNumber,Px, Py, Pz, Fx, Fy,Fz) 

nRows=nNodes 

nCols=7

Where all units are expressed in mm and the forces are left as variables

<img width="780" height="810" alt="image" src="https://github.com/user-attachments/assets/5cfb865f-253f-4f9c-a433-48bb3dc9edf8" />

The epsilon variables will be the one that are optimized by the system so are the degrees of freedom of the system. For example if a coordinate is defined as 100+epsilon1 when we start the optimization algorithm the system will try to find a better solution moving the coordinates around 100. We can also limit the value of the epsilon in order to avoid moving away from the desired position too much.

Not all coordinates will have a degree of freedom when optimizing the chassis because more variables we add and more arriving to convergence will be difficult for the system and more time will take (with 10 optimization parameter the system runs for 4/5 hours on a good pc). Furthermore we want to keep certain node coordinates where they are to respect certain rules, for the kinematic points…

F1 and F2 are the two forces applied to the suspensions and will be of value 1 to avoid numerical problems, with an higher value the program may not work. Knowing that deformation is linear in respect to the force applied if we want to evaluate the deformation of the chassis with 100N of force applied we just mutiply the obtained displacement value by 100.

F1 is a force used to apply the mohr maxwell method so during the simulation will be equal to 0. Depending on the position we place the F1 force we will calculate a different displacement. For example in the image above we are calculating the z displacement of node 25.

<img width="437" height="451" alt="image" src="https://github.com/user-attachments/assets/5560fd38-a843-47bf-acac-a7d817a3e92b" />

We create a vector containg the reaction forces that needs to be separated from the hyperstatic reaction forces. So the first vector will contain only the reaction forces that make the system isostatic. Currently the script works with when the system is once externally hyperstatic but if you need more than one it is a simple change of the code.

We then set some parameters. The density of the tubes used as an equivalence to the plates is set to 0 because their mass makes almost no difference while comparing different chassis.

### Beam matrix

Beam matrix is a table that contains information about all the tubes that composes the chassis. For example tube 7 is the tube that connects node number 4 and node number 8 and is a tube of 28x1.5mm caracterized as a beam. The second column indicates the starting point of the beam and the third column the ending point. Usually i have used the node with the lowest number as the start but the script works just fine reversing the start/end order.

We can have “Beam”, “Susp” and “Plates”.

“Susp” tubes can only transmit normal actions

“Plates” name is just used to color them in a different color in 3D plots.

<img width="411" height="643" alt="image" src="https://github.com/user-attachments/assets/c07045b3-a8b7-4471-bc8f-f91e66fa9f66" />

### 3D plot
In the 3D plot are shown some cyan arrows that indicates where the external forces and reaction forces are applied. This arrows do not represent the direction of the force but only the position.

<img width="1032" height="664" alt="image" src="https://github.com/user-attachments/assets/7454fb10-7a01-4fc8-b5aa-025c6973e648" />

### Functions to calculate beam parameters

This are some of the function that will be useful later in the code
<img width="1227" height="600" alt="image" src="https://github.com/user-attachments/assets/35cc6265-4ee6-4131-874c-27caa746cd42" />







