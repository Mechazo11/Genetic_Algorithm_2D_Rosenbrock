# Genetic_Algorithm_2D_Rosenbrock

Demonstrates usage of a Genetic Algorithm in solving a 2D unconstrained optimization problem. A standard Rosenbrock problem is chosen as an illustrative example

Version -- 3.5
Updated on -- 02/11/2023

Written in MATLAB Version -- 2019b

## TODO
1. find_mate_with_replacement and find_mate_without_replacement are hardcoded for 2 variable problems. This will be fixed in version 4 

Known dependencies
1. Statistics Toolbox -- randsample

## How to use
1. Run GA_rosenbrock.mlx
2. Choose N as a even number
3. Set "num_iter" to run Genetic Algorithm for "num_iter" generations
3. Supply U and L, the upper and lower bound values for x1,x2....xn
4. Choose number of bits to represent design variables, bit_count = [8,16,32]. I tested with 8 and 16
5. Input target_x1 = known global optimium location for variable x1
6. Input target_x2 = known global optimium location for variable x2
7. "Table" and "plot" control variables, takes a value between 0 and 1 [TODO need to check this again]

## Features
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
2nd year Ph.D. student, Dept. of MIE, Louisiana State University<br/>
Email: akamal4@lsu.edu
