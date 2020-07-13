This repository contains the Matlab implementations for the following multi-target filtering/tracking algorithms

- Folder PMBM contains the implementations of the Poisson multi-Bernoulli mixture (PMBM) filter and the multi-Bernoulli mixture (MBM) filter described in

A. F. García-Fernández, J. L. Williams, K. Granström, and L. Svensson, “Poisson multi-Bernoulli mixture filter: direct derivation and implementation,” IEEE Transactions on Aerospace and Electronic Systems, vol. 54, no. 4, pp. 1883–1901, Aug. 2018.

J. L. Williams, "Marginal multi-bernoulli filters: RFS derivation of MHT, JIPDA, and association-based member," in IEEE Transactions on Aerospace and Electronic Systems, vol. 51, no. 3, pp. 1664-1687, July 2015.

A. F. García-Fernández, Y. Xia , K. Granström, L. Svensson, J. L. Williams, "Gaussian implementation of the multi-Bernoulli mixture filter", in Proceedings of the 22nd International conference on Information Fusion, 2019.

In order to run the filters, execute PoissonMBMtarget_filter.m for the PMBM filter MBMtarget_filter.m for the MBM filter

- Folder TPHD contains the implementations of the trajectory probability hypothesis density (TPHD) filter and the trajectory cardinality PHD (TCPHD) filter for sets of trajectories in

A. F. García-Fernández and L. Svensson, “Trajectory PHD and CPHD filters”, IEEE Transactions on Signal Processing, vol. 67, no. 22, pp. 5702-5714,Nov. 2019.

Á. F. García-Fernández, L. Svensson and M. R. Morelande, "Multiple Target Tracking Based on Sets of Trajectories," in IEEE Transactions on Aerospace and Electronic Systems, vol. 56, no. 3, pp. 1685-1707, June 2020.

In order to run the filters, execute GM_TPHD_filter.m for the TPHD filter GM_TCPHD_filter.m for the TCPHD filter

- Folder CD MTT filters contains the implementations of the continuous-discrete PMBM, continuous-discrete PHD, and continuous-discrete CPHD filters described in

A. F. García-Fernández, S. Maskell, "Continuous-discrete multiple target filtering: PMBM, PHD and CPHD filter implementations," IEEE Transactions on Signal Processing, vol. 68, pp. 1300-1314, 2020.


- Evaluation of the multi-target filters is based on the generalised optimal subpattern-assignment (GOSPA) and its decomposition into localisation errors for properly detected targets, and costs for false and missed targets.

A. S. Rahmathullah, A. F. García-Fernández, and L. Svensson, “Generalized optimal sub-pattern assignment metric,” in 20th International Conference on Information Fusion, 2017.
A. F. García-Fernández, and L. Svensson, "Spooky effect in optimal OSPA estimation and how GOSPA solves it," in 22nd International Conference on Information Fusion, 2019.
Video on GOSPA: https://www.youtube.com/watch?v=M79GTTytvCM

- Evaluation of multi-target trackers (filters that estimate a set of trajectories) is based on the LP trajectory metric for sets of trajectories and its decomposition into localisation errors for properly detected targets, and costs for false, missed targets, and track switches.

Á. F. García-Fernández, A. S. Rahmathullah and L. Svensson, "A Metric on the Space of Finite Sets of Trajectories for Evaluation of Multi-Target Tracking Algorithms," in IEEE Transactions on Signal Processing, vol. 68, pp. 3917-3928, 2020.

- Open access versions of the above papers can be found in https://www.liverpool.ac.uk/electrical-engineering-and-electronics/staff/angel-garcia-fernandez/publications/

- A relevant online course on multiple target tracking is provided here:

https://www.youtube.com/channel/UCa2-fpj6AV8T6JK1uTRuFpw
