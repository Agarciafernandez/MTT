This repository contains the Matlab implementations of Bayesian multi-target filtering and tracking algorithms. The available algorithms and metrics are organised in the following folders.

# Performance evaluation

Performance evaluation of the multi-target filters, which estimate the set of targets at each time step, is based on the generalised optimal subpattern-assignment (GOSPA) metric and its decomposition into localisation errors for properly detected targets, and costs for false and missed targets  [[12]](https://ieeexplore.ieee.org/document/8009645)[[13]](https://ieeexplore.ieee.org/abstract/document/9011259)[[14]](https://www.youtube.com/watch?v=M79GTTytvCM). The GOSPA metric code is provided in folder "GOSPA code".

Performance evaluation of multi-target trackers (filters that estimate a set of trajectories) is based on the linear programming implementation of the trajectory GOSPA (T-GOSPA) metric (GOSPA metric for sets of trajectories) and its decomposition into localisation errors for properly detected targets, and costs for false, missed targets, and track switches [[15]](https://ieeexplore.ieee.org/document/9127194). The T-GOSPA metric code is provided in folder "Trajectory metric". This folder also contains the time-weighted T-GOSPA metric proposed in [[26]](https://ieeexplore.ieee.org/document/9626977).


# Folder "PMBM"

This folder contains the implementations of the Poisson multi-Bernoulli mixture (PMBM) filter [[1]](https://ieeexplore.ieee.org/document/7272821)[[2]](https://ieeexplore.ieee.org/document/8289337), the multi-Bernoulli mixture (MBM) filter [[3]](https://ieeexplore.ieee.org/abstract/document/9011346), and (track-oriented) Poisson multi-Bernoulli (PMB) [[1]](https://ieeexplore.ieee.org/document/7272821).

To run the PMBM filter, execute PMBMtarget_filter.m.
To run the PMB filter, execute PMBtarget_filter.m.
To run the MBM filter, execute MBMtarget_filter.m.

Sequential track estimation with the above filters can be achieved by linking target states estimates from the same Bernoulli component, which is uniquely identified by a start time and measurement.  This information can be made explicit in the posterior via auxiliary variables [[4]](https://ieeexplore.ieee.org/document/9169859). Note that Bayesian track formation is obtained by estimating the set of trajectories directly from the posterior density on the set of trajectories, not by sequentially linking target state estimates [[5]](https://ieeexplore.ieee.org/document/8731733).

To run the PMBM filter with sequential track formation, execute PMBMtarget_filter_tracks_all.m.
To run the PMB filter with sequential track formation, execute PMBtarget_filter_tracks_all.m.

# Folder "CD MTT filters"

This folder contains implementations of continuous-discrete multi-target filters. That is, these filters that have been designed for multi-target dynamics given in continuous time (including target births, single-target dynamics and target deaths) and measurements taken at known discrete time steps. These filters are specially suitable for non-uniform sampling times since the target dynamics should not be independent of how sensors take samples from the scene.

In particular, this folder contains the implementations of the continuous-discrete PMBM filter, the continuous-discrete probability hypothesis density (PHD) filter and continuous-discrete cardinality probability hypothesis density (CPHD) filter for the Wiener velocity model proposed in [[6]](https://ieeexplore.ieee.org/document/8964451):

To run the continuous-discrete PMBM filter (Wiener velocity), execute PMBMtarget_cd_filter.m.
To run the continuous-discrete PHD filter (Wiener velocity), execute GMPHD_cd_filter.m.
To run the continuous-discrete CPHD filter (Wiener velocity), execute GMCPHD_cd_filter.m.

This folder also contains the implementations of the continuous-discrete PMBM, continuous-discrete PMB, continuous-discrete PMB, continuous-discrete PHD, and continuous-discrete CPHD filters for general linear stochastic differential equations (SDEs) described in [[24]](https://ieeexplore.ieee.org/document/10857380):

To run the continuous-discrete PMBM filter (linear SDE), execute PMBMtarget_cd_filter_linear_SDE.m.
To run the continuous-discrete PMB filter (linear SDE), execute PMBtarget_cd_filter_linear_SDE.m.
To run the continuous-discrete PHD filter (linear SDE), execute GMPHD_cd_filter_linear_SDE.m.
To run the continuous-discrete CPHD filter (linear SDE), execute GMCPHD_cd_filter_linear_SDE.m.

This folder also contains the implementations of the continuous-discrete PMBM, continuous-discrete PMB, continuous-discrete PHD, and continuous-discrete CPHD filters for non-linear SDEs described in [[24]](https://ieeexplore.ieee.org/document/10857380):

To run the continuous-discrete PMBM filter (non-linear SDE), execute PMBMtarget_cd_filter_nonlinear_SDE_Taylor.m.
To run the continuous-discrete PMB filter (non-linear SDE), execute PMBtarget_cd_filter_nonlinear_SDE_Taylor.m.
To run the continuous-discrete PHD filter (non-linear SDE), execute GMPHD_cd_filter_nonlinear_SDE_Taylor.m.
To run the continuous-discrete CPHD filter (non-linear SDE), execute GMCPHD_cd_filter_nonlinear_SDE_Taylor.m.

# Folder "TPHD"

This folder contains the implementations of the trajectory probability hypothesis density (TPHD) filter and the trajectory cardinality PHD (TCPHD) filter for sets of trajectories in [[7]](https://ieeexplore.ieee.org/document/8846723).

To run the filters, execute GM_TPHD_filter.m and GM_TCPHD_filter.m

# Folder "TPMBM filter"

This folder contains the implementations of the trajectory PMBM (TPMBM) filter [[8]](https://ieeexplore.ieee.org/document/8455849)[[9]](https://ieeexplore.ieee.org/document/10799204), trajectory MBM (TMBM) filter [[10]](https://isif.org/media/multiscan-implementation-trajectory-poisson-multi-bernoulli-mixture-filter), trajectory PMB (TPMB) filter [4] and trajectory MB (TMB) filter [[11]](https://ieeexplore.ieee.org/document/9190554). Each of these filters can be run to estimate the set of alive trajectories or the set of all trajectories at each time step (running a different file).

To run the filters, execute TPMBM_all_filter.m, TPMB_all_filter.m, TMBM_all_filter.m, TMB_all_filter.m,  TPMBM_alive_filter.m, TPMB_alive_filter.m, TMBM_alive_filter.m, TMB_alive_filter.m.

# Folder "OOS TPMBM filter"

This folder contains the implementations of the continuous-discrete TPMBM and continuous-discrete  TPMB filters with out-of-sequence measurements [[16]](https://ieeexplore.ieee.org/document/9502575). That is, multi-target dynamics (including target births, single-target dynamics and target deaths) are modelled in continuous time, and we receive measurements at known discrete time steps. We can also receive measurements that happened before the last processed measurement, so this measurement is out-of-sequence.

To run the filters, execute TPMBM_cd_all_filter_oos.m and TPMB_cd_all_filter_oos.m.

# Folder "Tree PMBM - Spawning"

This folder contains the implementations of the Tree PMBM and Tree MBM filters for multiple target tracking with spawning [[17]](https://ieeexplore.ieee.org/document/9754270).

To run the filters, execute TrPMBM_all_filter.m and TrMBM_all_filter.m. A PMBM filter with spawning and sequential track formation can be run with PMBMtarget_filter_tracks_all_spawning.m.

# Folder "Non-linear MTT"

This folder contains the implementation of PMBM, PMB, TPMBM and TPMB filters for non-linear, non-Gaussian measurement models and non-constant probability of detection.
In particular, the filters perform the updates using the iterated posterior linearisation filter (IPLF) for conditional moments [[18]](https://ieeexplore.ieee.org/document/7153566)[[19]](https://ieeexplore.ieee.org/document/8260875) and with the improvement of normalising constant approximation [[20]](https://ieeexplore.ieee.org/document/9409712). The TPMBM and TPMB implementations are explained in [[21]](https://ieeexplore.ieee.org/document/9841396). The implementations make use of a range-bearings sensors with von-Mises Fisher-distributed measurements. The IPLF for this type of measurement is explained in [[22]](https://ieeexplore.ieee.org/document/8691412).

To run the filters, execute TPMBM_all_filter_range_bearing.m, TPMB_all_filter_range_bearing.m, TPMBM_alive_filter_range_bearing.m, TPMB_alive_filter_range_bearing.m, PMBM_filter_range_bearing.m , PMB_filter_range_bearing.m

# Folder "PMBM arbitrary clutter"

This folder contains the point-target PMBM and PMB implementations with negative binomial clutter in [[23]](https://ieeexplore.ieee.org/document/10130623).

To run the filters, execute "PMBMtarget_filter_nb_clutter" and "PMBtarget_filter_nb_clutter"

The code for the extended target example in [[23]](https://ieeexplore.ieee.org/document/10130623) can be found at

https://github.com/yuhsuansia/Extented-target-PMBM-filter-independent-clutter-sources

# Folder "Distributed MTT"

This folder contains the implementation of the distributed PMBM/PMB filters based on the generalised covariance intersection (GCI) fusion rule for PMB densities in [[25]](https://ieeexplore.ieee.org/document/11334185).

To run the filters, execute DPMB_GCI_filtering.m.


# Online course

A relevant online course on Bayesian multiple target tracking by Chalmers University is available in Youtube and edX:

https://www.youtube.com/channel/UCa2-fpj6AV8T6JK1uTRuFpw

https://www.edx.org/course/multi-object-tracking-for-automotive-systems

# References

[[1] J. L. Williams, "Marginal multi-Bernoulli filters: RFS derivation of MHT, JIPDA, and association-based MeMBer," in IEEE Transactions on Aerospace and Electronic Systems, vol. 51, no. 3, pp. 1664-1687, July 2015.](https://ieeexplore.ieee.org/document/7272821)

[[2] A. F. García-Fernández, J. L. Williams, K. Granström, and L. Svensson,  Poisson multi-Bernoulli mixture filter: direct derivation and implementation,  IEEE Transactions on Aerospace and Electronic Systems, vol. 54, no. 4, pp. 1883 1901, Aug. 2018.](https://ieeexplore.ieee.org/document/8289337)

[[3] A. F. García-Fernández, Y. Xia , K. Granström, L. Svensson, J. L. Williams, "Gaussian implementation of the multi-Bernoulli mixture filter", in Proceedings of the 22nd International conference on Information Fusion, 2019.](https://ieeexplore.ieee.org/abstract/document/9011346)

[[4] A. F. García-Fernández, L. Svensson, J. L. Williams, Y. Xia and K. Granström, "Trajectory Poisson Multi-Bernoulli Filters," in IEEE Transactions on Signal Processing, vol. 68, pp. 4933-4945, 2020.](https://ieeexplore.ieee.org/document/9169859)

[[5] A. F. García-Fernández, L. Svensson and M. R. Morelande, "Multiple Target Tracking Based on Sets of Trajectories," in IEEE Transactions on Aerospace and Electronic Systems, vol. 56, no. 3, pp. 1685-1707, June 2020.](https://ieeexplore.ieee.org/document/8731733)

[[6] A. F. García-Fernández, S. Maskell, "Continuous-discrete multiple target filtering: PMBM, PHD and CPHD filter implementations," IEEE Transactions on Signal Processing, vol. 68, pp. 1300-1314, 2020.](https://ieeexplore.ieee.org/document/8964451)

[[7] A. F. García-Fernández and L. Svensson,  Trajectory PHD and CPHD filters , IEEE Transactions on Signal Processing, vol. 67, no. 22, pp. 5702-5714, Nov. 2019.](https://ieeexplore.ieee.org/document/8846723)

[[8] K. Granström, L. Svensson, Y. Xia, J. Williams and A . F. García-Fernández, "Poisson Multi-Bernoulli Mixture Trackers: Continuity Through Random Finite Sets of Trajectories," 2018 21st International Conference on Information Fusion (FUSION), Cambridge, 2018.](https://ieeexplore.ieee.org/document/8455849)

[[9] K. Granström, L. Svensson, Y. Xia, J. Williams and A . F. García-Fernández, "Poisson Multi-Bernoulli Mixtures for Sets of Trajectories," in IEEE Transactions on Aerospace and Electronic Systems,vol. 61, no. 2, pp. 5178-5194, April 2025.](https://ieeexplore.ieee.org/document/10799204)

[[10] Y. Xia, K. Granström, L. Svensson, A. F. García-Fernández, and J. L. Wlliams,  "Multi-scan implementation of the trajectory Poisson multi-Bernoulli mixture filter",  Journal of Advances in Information Fusion. Special Issue on Multiple Hypothesis Tracking., vol. 14, no. 2, pp. 213 235, Dec. 2019.](https://isif.org/media/multiscan-implementation-trajectory-poisson-multi-bernoulli-mixture-filter)

[[11] A. F. García-Fernández, L. Svensson, J. L. Williams, Y. Xia, K. Granström,  "Trajectory multi-Bernoulli filters for multi-target tracking based on sets of trajectories,"  in 23rd International Conference on Information Fusion, 2020.](https://ieeexplore.ieee.org/document/9190554)

[[12] A. S. Rahmathullah, A. F. García-Fernández, and L. Svensson,  "Generalized optimal sub-pattern assignment metric,"  in 20th International Conference on Information Fusion, 2017.](https://ieeexplore.ieee.org/document/8009645)

[[13] A. F. García-Fernández, and L. Svensson, "Spooky effect in optimal OSPA estimation and how GOSPA solves it," in 22nd International Conference on Information Fusion, 2019.](https://ieeexplore.ieee.org/abstract/document/9011259)

[14] L. Svensson, Presentation on GOSPA metric: https://www.youtube.com/watch?v=M79GTTytvCM

[[15] A. F. García-Fernández, A. S. Rahmathullah and L. Svensson, "A Metric on the Space of Finite Sets of Trajectories for Evaluation of Multi-Target Tracking Algorithms," in IEEE Transactions on Signal Processing, vol. 68, pp. 3917-3928, 2020.](https://ieeexplore.ieee.org/document/9127194)

[[16] A. F. García-Fernández and W. Yi, "Continuous-Discrete Multiple Target Tracking With Out-of-Sequence Measurements," in IEEE Transactions on Signal Processing, vol. 69, pp. 4699-4709, 2021](https://ieeexplore.ieee.org/document/9502575)

[[17] A. F. García-Fernández and L. Svensson, "Tracking Multiple Spawning Targets Using Poisson Multi-Bernoulli Mixtures on Sets of Tree Trajectories," in IEEE Transactions on Signal Processing, vol. 70, pp. 1987-1999, 2022.](https://ieeexplore.ieee.org/document/9754270)

[[18] A. F. García-Fernández, L. Svensson, M. R. Morelande and S. Särkkä, "Posterior Linearization Filter: Principles and Implementation Using Sigma Points," in IEEE Transactions on Signal Processing, vol. 63, no. 20, pp. 5561-5573, Oct.15, 2015.](https://ieeexplore.ieee.org/document/7153566)

[[19] F. Tronarp,  A. F. García-Fernández and S. Särkkä, "Iterative Filtering and Smoothing in Nonlinear and Non-Gaussian Systems Using Conditional Moments,"  in IEEE Signal Processing Letters, vol. 25, no. 3, pp. 408-412, March 2018.](https://ieeexplore.ieee.org/document/8260875)

[[20] A. F. García-Fernández,J. Ralph, P. Horridge, S. Maskell, "A Gaussian filtering method for multi-target tracking with nonlinear/non-Gaussian measurements" IEEE Transactions on Aerospace and Electronic Systems, vol. 57, no. 5, pp. 3539-3548, Oct. 2021.](https://ieeexplore.ieee.org/document/9409712)

[[21] A. F. García-Fernández, J. Ralph, P. Horridge, S. Maskell, "Gaussian trajectory PMBM filter with nonlinear measurements based on posterior linearisation" in Proceedings of the 25th International Conference on Information Fusion, 2022.](https://ieeexplore.ieee.org/document/9841396)

[[22]  A. F. García-Fernández, F. Tronarp and S. Särkkä, "Gaussian Target Tracking With Direction-of-Arrival von Mises Fisher Measurements," in IEEE Transactions on Signal Processing, vol. 67, no. 11, pp. 2960-2972, June, 2019.](https://ieeexplore.ieee.org/document/8691412)

[[23] A. F. García-Fernández, Y. Xia, L. Svensson, "Poisson multi-Bernoulli mixture filter with general target-generated measurements and arbitrary clutter", IEEE Transactions on Signal Processing, vol. 71, pp. 1895-1906, 2023.](https://ieeexplore.ieee.org/document/10130623)

[[24] A. F. García-Fernández, S. Särkkä, "Gaussian multi-target filtering with target dynamics driven by a stochastic differential equation", in IEEE Transactions on Signal Processing, vol. 73, pp. 664-675, 2025.](https://ieeexplore.ieee.org/document/10857380)

[[25] Á. F. García-Fernández and G. Battistelli, "Distributed Poisson Multi-Bernoulli Filtering via Generalized Covariance Intersection," in IEEE Transactions on Signal Processing, vol. 74, pp. 246-257, 2026.](https://ieeexplore.ieee.org/document/11334185)

[[26] Á. F. García-Fernández, A. S. Rahmathullah and L. Svensson, "A time-weighted metric for sets of trajectories to assess multi-object tracking algorithms" in Proceedings of the 24th International Conference on Information Fusion, 2021.](https://ieeexplore.ieee.org/document/9626977)
