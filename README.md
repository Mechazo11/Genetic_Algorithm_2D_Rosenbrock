# Genetic_Algorithm_2D_Rosenbrock

Genetic Algorithm that solves 2D Rosenbrock problem. Look into GA_rosenbrock.mlx for details.
Written in MATLAB Version -- 2019b

Known dependencies
1. Statistics Toolbox -- randsample

Usage
1. Run GA_rosenbrock.mlx
2. Values of candidate, N must be even number.

Features
* Crossover probabilty randomized -- (To be implemented)
* Check N as even -- (To be implemented) 
* Mutation enabled
* Mutation occurs at child level.
* Bit-Flip Mutation operator - modified approach, look at function bit_flip(child_array)
* Random Building Block Operator -- (To be implemented) 
* Elitisim (Copy fittest individual to next generation)
* Survival of the fittest (drop a candidate with poorest fit to bring N+1 candidate pool to N candidate pool)
* Allowed gene length = 8, 16 (can go up as power of 2 but only tested 8 and 16)
* Starting population randomly created
* Auto termination -- (To be implemented)

For questions, suggestions please feel free to contact me

With best,
Azmyin Md. Kamal,
Graduate Research Assistant,
M.Sc. in Mechanical Engineering,
Univ. of Louisiana at Lafayette,
Email: c00441440@louisiana.edu
