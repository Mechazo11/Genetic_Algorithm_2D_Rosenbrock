# Genetic_Algorithm_2D_Rosenbrock

Genetic Algorithm that solves 2D Rosenbrock problem. Look into GA_rosenbrock.mlx for details.

Version -- 3.0
Updated on -- 11/08/2020

Written in MATLAB Version -- 2019b

Known dependencies
1. Statistics Toolbox -- randsample

Usage
1. Run GA_rosenbrock.mlx
2. Choose N as a even number
3. Choose num_iter to run Genetic Algorithm for "num_iter" generations
3. Supply U and L, the upper and lower bound values for x1,x2....xn
4. Choose number of bits to represent design variables, bit_count = [8,16,32]. I tested with 8 and 16
5. Input target_x1 = known global optimium location for variable x1
6. Input target_x2 = known global optimium location for variable x2
7. Table and plot control variables can take a value between 0 and 1


Features
* Single point cross over
* Bit flip mutation operator
* Mate pool -- choose with replacement (i.e. a candiate can be repeated)
* Elitisim (copy strongest candiate to next generation)
* Survival of the fittest, Gene enconding = 8bits, 16bit
* Starting population randomly created
* Script can be adopted for any two variable design problem

For questions, suggestions please feel free to contact me

With best regards,
Azmyin Md. Kamal,<br/>
Graduate Research Assistant,<br/>
M.Sc. in Mechanical Engineering,<br/>
The University of Louisiana at Lafayette,<br/>
Email: c00441440@louisiana.edu
