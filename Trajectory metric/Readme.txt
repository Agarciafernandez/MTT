This folder contains the implementations of two trajectory generalised optimal subpattern assignment (GOSPA) metrics (metrics for sets of trajectories) to measure the accuracy of multi-target tracking algorithms. 
These metrics are principled mathematical metrics, which meet the identity, symmetry and triangle inequality properties. 
They have been designed to penalise localisation errors for properly detected targets, missed targets, false targets and track switches. 
The decomposition of the trajectory GOSPA error in terms of the above components is also provided in the code below. This decomposition is important to understand and analyse the performance of the evaluated multi-target tracking algorithms. 

- examplesForTrajectoryMetric.m runs a demo with the trajectory GOSPA metric in [1] based on linear programming.
- LPTrajMetric_cluster.m contains the trajectory GOSPA implementation of [1], based on linear programming and clustering.
- examplesForTWTrajectoryMetric runs a demo with the time-weighted trajectory GOSPA metric in [2], also based on linear programming.
- TimeWeightedLPTrajMetric_cluster.m contains the implementation of the time-weighted trajectory GOSPA metric [2], based on linear programming and clustering.

[1] Á. F. García-Fernández, A. S. Rahmathullah and L. Svensson,"A Metric on the Space of Finite Sets of Trajectories for Evaluation of 
Multi-Target Tracking Algorithms," in IEEE Transactions on Signal Processing, vol. 68, pp. 3917-3928, 2020.

[2] Á. F. García-Fernández, A. S. Rahmathullah and L. Svensson, "A time-weighted metric for sets of trajectories to assess multi-object
tracking algorithms" in Proceedings of the 24th International Conference on Information Fusion, 2021.