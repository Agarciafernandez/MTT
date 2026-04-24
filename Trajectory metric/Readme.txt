This folder contains the Matlab implementations of the T-GOSPA metric [1] and the time-weighted T-GOSPA metric [2] for sets of trajectories to measure the accuracy of multi-target tracking algorithms. These metrics penalise the localisation error for properly detected targets, the number of missed targets, the number of false targets, and the number of track switches.

- examplesForTrajectoryMetric.m runs a demo with the T-GOSPA metric [1].
- LPTrajMetric_cluster.m contains the T-GOSPA metric implementation of [1], based on linear programming and clustering.
- examplesForTWTrajectoryMetric runs a demo with the time-weighted T-GOSPA metric [2].
- TimeWeightedLPTrajMetric_cluster.m contains the time-weighted T-GOSPA metric implementation of [2], based on linear programming and clustering.

The folder also contains the Matlab implementation of the T-GOSPA quasi-metric (q-metric) [3]. The T-GOSPA q-metric follows the same principles as the T-GOSPA metric but enables different costs for missed and false targets.

- examplesForTrajectoryQMetric.m runs a demo with the T-GOSPA q-metric [3].

- LPTrajQMetric_cluster.m contains the T-GOSPA q-metric implementation of [3], based on linear programming and clustering.


[1] Á. F. García-Fernández, A. S. Rahmathullah and L. Svensson,"A Metric on the Space of Finite Sets of Trajectories for Evaluation of 
Multi-Target Tracking Algorithms," in IEEE Transactions on Signal Processing, vol. 68, pp. 3917-3928, 2020.

[2] Á. F. García-Fernández, A. S. Rahmathullah and L. Svensson, "A time-weighted metric for sets of trajectories to assess multi-object
tracking algorithms" in Proceedings of the 24th International Conference on Information Fusion, 2021.

[3] A. F. García-Fernández, J. Gu, L. Svensson, Y. Xia, J. Krejcí, O. Kost, O. Straka “GOSPA and T-GOSPA quasi-metrics for evaluation of multi-object tracking algorithms,” accepted in IEEE Transactions on Aerospace and Electronic Systems, 2026.