This repository contains the Matlab implementations for the following multi-target filtering/tracking algorithms

- Folder PMBM contains the implementations of the Poisson multi-Bernoulli mixture (PMBM) filter and the multi-Bernoulli mixture (MBM) filter described in

A. F. García-Fernández, J. L. Williams, K. Granström, and L. Svensson, “Poisson multi-Bernoulli mixture filter: direct derivation and implementation,” IEEE Transactions on Aerospace and Electronic Systems, vol. 54, no. 4, pp. 1883–1901, Aug. 2018.
A. F. García-Fernández,Y. Xia , K. Granström, L. Svensson, J. L. Williams, "Gaussian implementation of the multi-Bernoulli mixture filter", in Proceedings of the 22nd International conference on Information Fusion, 2019.

In order to run the filters, execute
PoissonMBMtarget_filter.m for the PMBM filter
MBMtarget_filter.m for the MBM filter

- Folder TPHD contains the implementations of the trajectory probability hypothesis density (TPHD) filter and the trajectory cardinality PHD (TCPHD) filter for sets of trajectories in

A. F. García-Fernández and L. Svensson, “Trajectory PHD and CPHD filters”, IEEE Transactions on Signal Processing, vol. 67, no. 22, pp. 5702-5714,Nov. 2019.
A. F. García-Fernández, L. Svensson, and M. R. Morelande, “Multiple target tracking based on sets of trajectories,” IEEE Transactions on Aerospace and Electronic Systems. https://arxiv.org/abs/1605.08163

In order to run the filters, execute
GM_TPHD_filter.m for the TPHD filter
GM_TCPHD_filter.m for the TCPHD filter


- Evaluation of the multi-target filters is based on the generalised optimal subpattern-assignment (GOSPA) and its decomposition into localisation errors for properly detected targets, and costs for false and missed targets.

A. S. Rahmathullah, A. F. García-Fernández, and L. Svensson, “Generalized optimal sub-pattern assignment metric,” in 20th International Conference on Information Fusion, 2017.
A. F. García-Fernández, and L. Svensson, "Spooky effect in optimal OSPA estimation and how GOSPA solves it," in 22nd International Conference on Information Fusion, 2019.
Video on GOSPA: https://www.youtube.com/watch?v=M79GTTytvCM

- Evaluation of multi-target trackers is based on the LP trajectory metric for sets of trajectories and its decomposition into localisation errors for properly detected targets, and costs for false, missed targets, and track switches.

A. S. Rahmathullah, A. F. García-Fernández, and L. Svensson, “A metric on the space of ?nite sets of trajectories for evaluation of multi-target tracking algorithms,” 2016. http://arxiv.org/abs/1605.01177

- A relevant online course on multiple target tracking is provided here:

https://www.youtube.com/channel/UCa2-fpj6AV8T6JK1uTRuFpw
