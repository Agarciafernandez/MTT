This repository contains the Matlab implementations for the following multi-target filtering/tracking algorithms:

- Folder PMBM contains the implementations of the Poisson multi-Bernoulli mixture (PMBM) filter and the multi-Bernoulli mixture (MBM) filter described in [1]-[3].


In order to run the filters, execute PMBMtarget_filter.m for the PMBM filter MBMtarget_filter.m for the MBM filter

PMBMtarget_filter_tracks_all.m runs the PMBM filter with sequential track formation, linking target states estimates from the same Bernoulli component, which is uniquely identified by a start time and measurement. This information can be made explicit in the posterior via auxiliary variables [4]. Note that Bayesian track formation is obtained via densities on sets of trajectories, not linking target state estimates [5].

- Folder CD MTT filters contains the implementations of the continuous-discrete PMBM, continuous-discrete PHD, and continuous-discrete CPHD filters described in [6].

- Folder TPHD contains the implementations of the trajectory probability hypothesis density (TPHD) filter and the trajectory cardinality PHD (TCPHD) filter for sets of trajectories in [7].

In order to run the filters, execute GM_TPHD_filter.m and GM_TCPHD_filter.m

- Folder TPMBM filter contains the implementations of the trajectory PMBM (TPMBM) filter [8][9], trajectory MBM (TMBM) filter [10], trajectory PMB (TPMB) filter [4] and trajectory MB (TMB) filter [11]. Each of these filters can be run to estimate the set of alive trajectories or the set of all trajectories at each time step (running a different file).



- Evaluation of the multi-target filters is based on the generalised optimal subpattern-assignment (GOSPA) and its decomposition into localisation errors for properly detected targets, and costs for false and missed targets  [12][13][14].


- Evaluation of multi-target trackers (filters that estimate a set of trajectories) is based on the LP trajectory metric for sets of trajectories and its decomposition into localisation errors for properly detected targets, and costs for false, missed targets, and track switches [15].


- Open access versions of the above papers can be found in https://www.liverpool.ac.uk/electrical-engineering-and-electronics/staff/angel-garcia-fernandez/publications/

- A relevant online course on multiple target tracking is provided here:

https://www.youtube.com/channel/UCa2-fpj6AV8T6JK1uTRuFpw

REFERENCES

[1] A. F. García-Fernández, J. L. Williams, K. Granström, and L. Svensson, “Poisson multi-Bernoulli mixture filter: direct derivation and implementation,” IEEE Transactions on Aerospace and Electronic Systems, vol. 54, no. 4, pp. 1883–1901, Aug. 2018.

[2] J. L. Williams, "Marginal multi-bernoulli filters: RFS derivation of MHT, JIPDA, and association-based member," in IEEE Transactions on Aerospace and Electronic Systems, vol. 51, no. 3, pp. 1664-1687, July 2015.

[3] A. F. García-Fernández, Y. Xia , K. Granström, L. Svensson, J. L. Williams, "Gaussian implementation of the multi-Bernoulli mixture filter", in Proceedings of the 22nd International conference on Information Fusion, 2019.

[4] Á. F. García-Fernández, L. Svensson, J. L. Williams, Y. Xia and K. Granström, "Trajectory Poisson Multi-Bernoulli Filters," in IEEE Transactions on Signal Processing, vol. 68, pp. 4933-4945, 2020.

[5] Á. F. García-Fernández, L. Svensson and M. R. Morelande, "Multiple Target Tracking Based on Sets of Trajectories," in IEEE Transactions on Aerospace and Electronic Systems, vol. 56, no. 3, pp. 1685-1707, June 2020.

[6] A. F. García-Fernández, S. Maskell, "Continuous-discrete multiple target filtering: PMBM, PHD and CPHD filter implementations," IEEE Transactions on Signal Processing, vol. 68, pp. 1300-1314, 2020.

[7] A. F. García-Fernández and L. Svensson, “Trajectory PHD and CPHD filters”, IEEE Transactions on Signal Processing, vol. 67, no. 22, pp. 5702-5714,Nov. 2019.

[8] K. Granström, L. Svensson, Y. Xia, J. Williams and Á. F. García-Fernández, "Poisson Multi-Bernoulli Mixture Trackers: Continuity Through Random Finite Sets of Trajectories," 2018 21st International Conference on Information Fusion (FUSION), Cambridge, 2018.

[9] K. Granström, L. Svensson, Y. Xia, J. Williams and Á. F. García-Fernández, "Poisson Multi-Bernoulli Mixtures for Sets of Trajectories," https://arxiv.org/abs/1912.08718

[10] Y. Xia, K. Granström, L. Svensson, A. F. García-Fernández, and J. L. Wlliams, “Multi-scan implementation of the trajectory Poisson multi-Bernoulli mixture filter,” Journal of Advances in Information Fusion. Special Issue on Multiple Hypothesis Tracking., vol. 14, no. 2, pp. 213–235, Dec. 2019.

[11] A. F. García-Fernández, L. Svensson, J. L. Williams, Y. Xia, K. Granström,  “Trajectory multi-Bernoulli filters for multi-target tracking based on sets of trajectories” in 23rd International Conference on Information Fusion, 2020.

[12] A. S. Rahmathullah, A. F. García-Fernández, and L. Svensson, “Generalized optimal sub-pattern assignment metric,” in 20th International Conference on Information Fusion, 2017.

[13] A. F. García-Fernández, and L. Svensson, "Spooky effect in optimal OSPA estimation and how GOSPA solves it," in 22nd International Conference on Information Fusion, 2019.

[14] L. Svensson, Presentation on GOSPA: https://www.youtube.com/watch?v=M79GTTytvCM

[15] Á. F. García-Fernández, A. S. Rahmathullah and L. Svensson, "A Metric on the Space of Finite Sets of Trajectories for Evaluation of Multi-Target Tracking Algorithms," in IEEE Transactions on Signal Processing, vol. 68, pp. 3917-3928, 2020.






