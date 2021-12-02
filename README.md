# Title
Wireless Sensor Network Deployment Optimization Based on Two Flower Pollination Algorithms
# DOI
10.1109/ACCESS.2019.2959949
# Abstract
For the wireless sensor networks (WSNs) heterogeneous node deployment optimization problem with obstacles in the monitoring area, two new flower pollination algorithms (FPA) are proposed to deploy the network. Firstly, an improved flower pollination algorithm (IFPA) is proposed based on FPA, aiming at the shortcomings of the convergence speed is slow and the precision is not high enough of FPA. The nonlinear convergence factor is designed to correct the scaling factor of FPA, the Tent chaotic map effectively maintains the diversity of the population in the late iteration, and a greedy crossover strategy is designed to assist the remaining individual search with better individuals. Secondly, based on FPA, a non-dominated sorting multi-objective flower pollination algorithm (NSMOFPA) is proposed. The external archive strategy and leader strategy are introduced, to solve the global pollination problem. The proposed crowding degree method and the introduced elite strategy effectively maintain the diversity of the population. Then, IFPA is applied to WSN deployment aiming at optimizing coverage rate, simulation experiments show that IFPA can obtain a higher coverage rate with shorter iterations, which can save network deployment costs. Finally, applying NSMOFPA to the WSN deployment with optimization objectives for coverage rate, node radiation overflow rate and energy consumption rate. The experimental results verify that NSMOFPA has a good optimization effect and can provide a better solution for WSN deployment.
# Environment
Matlab 2021a , If the code runs in error, install the symbolic math toolbox
# Running
init_data.m, IFPA_WSN.m, NSMOFPA.m
# Dataset
No
# Citation
If you use this codebase or any part of it for a publication, please cite:
Wang Z, Xie H, He D, et al. Wireless sensor network deployment optimization based on two flower pollination algorithms[J]. IEEE Access, 2019, 7: 180590-180608.
# Contact
fhzm1995@163.com
