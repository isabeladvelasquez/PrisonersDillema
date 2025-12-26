My name is Isabela Velasquez and here is where I'll keep the evolution of my project. 
\
Based on the Prisoner's Dillema through agents simulations in C, I am focused in understading the conditions which will lead to cooperation endurance. At first, I am building simulations that represent system without movement, such as total density of the network or just stopped agents. These nets are all non-deterministic, so they do not have the same density (agent per space) neither the same distribuition of strategies (C or D). 
Currently, I am adding movement to the dynamics and reproducing some results preexistents.
\
I began writing every code in brazilian portuguese, yet I already started to switch to english.
\
For now, below there is the description of the codes I uploaded.
\
\
L100_M0_RHO1_S42.c : This is the code I began with. It has a net with full density and time(NULL) as the seed. It was stablished a net with 10.000 spaces, all of them with agents which have their strategy decided with a 50% chance of being a cooperator (C) or a non cooperator (D). It can also make a image of the net in colors, with C as blue and D as red, at the beggining and at the end of the simulation, as well as print the density of cooperators in this same instants.
\
L100_M0_RHOvar.c : In this code you will find a simullation that covers 20 values of density, from 0.05 to 1, plus 0. In addition, I have created funtions to save the density of the net, the density of cooperators by the density of the net, as well as the evolution of the number of cooperators during the whole simulation. Also, you can print the net (with the active sites) just by adding the corresponding funtion at any time. All of this data will be saved on the folder that your code is running and you can allow ou deny the data saving by just changing the responsible defined constants at the beggining of the code.
\
repeat.sh : It's a code in bash that I use to repeat a code a "reps" number of times. It only needs the name of the code and the folder which the outputs will go, but you can change the number of repetitions to your liking.
\
analise.py : This code was built to analyze all the data I gathered using "repeat.sh" on L100_M0_RHOvar.c, which clarifies why this was written on python3. I decided to make one code for everything to keep my folders neat, so this single code has a funtion to each type of data I am analyzing. You may find as well some funtions to fit curves, which were used to make it clear the behavior of the population in some graphics. The libraries I used were Numpy, SciPy, Matplotlib.


