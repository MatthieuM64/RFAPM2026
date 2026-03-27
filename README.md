# Active Potts model under an external field

Codes used in the scientific publication: M. Karmakar, M. Mangeat, S. Chatterjee, H. Rieger, and R. Paul, <i>Field‐controlled interfacial transport and pinning in an active spin system</i>, submitted (2026). A preprint is available on arXiv.</br></br>

For each considered field (unidirectional, bidirection and random orientational), a C++ code to compute the numerical simulations of the microscopic model and a C++ code to compute the numerical solutions of the hydrodynamic equations are available in this repository. </br></br>
<b>Exportations:</b> density snapshots and profiles shown in the different figures of the paper.</br>
<b>Compile:</b> g++ filename.cpp -fopenmp -lgsl -lgslcblas -lm -O3 -s -o filename.out.</br>
<b>Run:</b> ./filename.out -parameter=value.</br>
<b>List of parameters for the numerical simulations</b>: beta, h, D, epsilon, rho0, LX, LY, init, tmax, ran, threads (details as comments in the code).</br>
<b>List of parameters for the hydrodynamic solutions</b>: beta, epsilon, rho0, h, LX, LY, init, dt, tmax, dx, ran, threads (details as comments in the code).</br>

## Active Potts model under a homogeneous unidirectional field

Minimal model with a homogeneous field in the right direction.</br>
<b>Files:</b> RFAPM_micro_unidirectional_omp.cpp; RFAPM_hydro_unidirectional_omp.cpp</br>

## Active Potts model under a bidirectional field

Model with a bidirectional field: right-oriented field in the left part, left-oriented field in the right part.</br>
<b>Files:</b> RFAPM_micro_bidirectional_omp.cpp; RFAPM_hydro_bidirectional_omp.cpp</br>

## Active Potts model under a random orientational field

Model with a random orientational field: the field orientation is random among the 4-states at each lattice site.</br>
<b>Files:</b> RFAPM_micro_random_omp.cpp; RFAPM_hydro_random_omp.cpp</br>
