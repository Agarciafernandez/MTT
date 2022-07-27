This repository contains the Matlab implementations for the following multi-target filtering/tracking algorithms:

- Folder "PMBM" contains the implementations of the Poisson multi-Bernoulli mixture (PMBM) filter [1][2], the multi-Bernoulli mixture (MBM) filter [3], and (track-oriented) Poisson multi-Bernoulli (PMB) [1].


To run the filters, execute PMBMtarget_filter.m for the PMBM filter MBMtarget_filter.m for the MBM filter

PMBMtarget_filter_tracks_all.m runs the PMBM filter with sequential track formation, linking target states estimates from the same Bernoulli component, which is uniquely identified by a start time and measurement. This information can be made explicit in the posterior via auxiliary variables [4]. Note that Bayesian track formation is obtained via densities on sets of trajectories, not linking target state estimates [5].

- Folder "CD MTT filters" contains the implementations of the continuous-discrete PMBM, continuous-discrete PHD, and continuous-discrete CPHD filters described in [6].

To run the filters, execute PoissonMBMtarget_cd_filter.m, GMPHD_cd_filter.m, GMCPHD_cd_filter.m.

- Folder "TPHD" contains the implementations of the trajectory probability hypothesis density (TPHD) filter and the trajectory cardinality PHD (TCPHD) filter for sets of trajectories in [7].

To run the filters, execute GM_TPHD_filter.m and GM_TCPHD_filter.m

- Folder "TPMBM filter" contains the implementations of the trajectory PMBM (TPMBM) filter [8][9], trajectory MBM (TMBM) filter [10], trajectory PMB (TPMB) filter [4] and trajectory MB (TMB) filter [11]. Each of these filters can be run to estimate the set of alive trajectories or the set of all trajectories at each time step (running a different file).

- To run the filters, execute TPMBM_all_filter.m, TPMB_all_filter.m, TMBM_all_filter.m, TMB_all_filter.m,  TPMBM_alive_filter.m, TPMB_alive_filter.m, TMBM_alive_filter.m, TMB_alive_filter.m. 

- Folder "OOS TPMBM filter" contains the implementations of the continuous-discrete TPMBM and continuous-discrete  TPMB filters with out-of-sequence measurements [16].

To run the filters, execute TPMBM_cd_all_filter_oos.m and TPMB_cd_all_filter_oos.m. 

- Folder "Tree PMBM - Spawning" contains the implementations of the Tree PMBM and Tree MBM filters for multiple target tracking with spawning [17].

To run the filters, execute TrPMBM_all_filter.m and TrMBM_all_filter.m. A PMBM filter with spawning and sequential track formation can be run with PMBMtarget_filter_tracks_all_spawning.m. 

- Folder "Non-linear MTT" contains the implementation of PMBM, PMB, TPMBM and TPMB filters for non-linear, non-Gaussian measurement models and non-constant probability of detection. 
In particular, the filters perform the updates using the iterated posterior linearisation filter (IPLF) for conditional moments [18][19] and improvement of normalising constant approximation [20]. The TPMBM and TPMB implementations are explained in [21]. The implementations make use of a range-bearings sensors with von-Mises Fisher-distributed measurements. The IPLF for this type of measurement is explained in [22].

To run the filters, execute TPMBM_all_filter_range_bearing.m, TPMB_all_filter_range_bearing.m, TPMBM_alive_filter_range_bearing.m, TPMB_alive_filter_range_bearing.m, PMBM_filter_range_bearing.m , PMB_filter_range_bearing.m 


- Evaluation of the multi-target filters is based on the generalised optimal subpattern-assignment (GOSPA) and its decomposition into localisation errors for properly detected targets, and costs for false and missed targets  [12][13][14].


- Evaluation of multi-target trackers (filters that estimate a set of trajectories) is based on the LP trajectory metric for sets of trajectories and its decomposition into localisation errors for properly detected targets, and costs for false, missed targets, and track switches [15].


- Open access versions of the above papers can be found in https://www.liverpool.ac.uk/electrical-engineering-and-electronics/staff/angel-garcia-fernandez/publications/

- A relevant online course on multiple target tracking by Chalmers University is provided in Youtube and edX:

https://www.youtube.com/channel/UCa2-fpj6AV8T6JK1uTRuFpw
https://www.edx.org/course/multi-object-tracking-for-automotive-systems 



REFERENCES

[1] J. L. Williams, "Marginal multi-bernoulli filters: RFS derivation of MHT, JIPDA, and association-based member," in IEEE Transactions on Aerospace and Electronic Systems, vol. 51, no. 3, pp. 1664-1687, July 2015.

[2] A. F. García-Fernández, J. L. Williams, K. Granström, and L. Svensson, “Poisson multi-Bernoulli mixture filter: direct derivation and implementation,” IEEE Transactions on Aerospace and Electronic Systems, vol. 54, no. 4, pp. 1883–1901, Aug. 2018.

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

[16] Á. F. García-Fernández and W. Yi, "Continuous-Discrete Multiple Target Tracking With Out-of-Sequence Measurements," in IEEE Transactions on Signal Processing, vol. 69, pp. 4699-4709, 2021

[17] Á. F. García-Fernández and L. Svensson, "Tracking Multiple Spawning Targets Using Poisson Multi-Bernoulli Mixtures on Sets of Tree Trajectories," in IEEE Transactions on Signal Processing, vol. 70, pp. 1987-1999, 2022.

[18] Á. F. García-Fernández, L. Svensson, M. R. Morelande and S. Särkkä, "Posterior Linearization Filter: Principles and Implementation Using Sigma Points," in IEEE Transactions on Signal Processing, vol. 63, no. 20, pp. 5561-5573, Oct.15, 2015.

[19] F. Tronarp, Á. F. García-Fernández and S. Särkkä, "Iterative Filtering and Smoothing in Nonlinear and Non-Gaussian Systems Using Conditional Moments,"  in IEEE Signal Processing Letters, vol. 25, no. 3, pp. 408-412, March 2018.

[20] A. F García-Fernández,J. Ralph, P. Horridge, S. Maskell, "A Gaussian filtering method for multi-target tracking with nonlinear/non-Gaussian measurements" IEEE Transactions on Aerospace and Electronic Systems, 2021, 57, 3539-3548.

[21] Á. F. García-Fernández, J. Ralph, P. Horridge, S. Maskell, "Gaussian trajectory PMBM filter with nonlinear measurements based on posterior linearisation" in Proceedings of the 25th International Conference on Information Fusion, 2022.

[22] Á. F. García-Fernández, F. Tronarp and S. Särkkä, "Gaussian Target Tracking With Direction-of-Arrival von Mises–Fisher Measurements," in IEEE Transactions on Signal Processing, vol. 67, no. 11, pp. 2960-2972, June, 2019.