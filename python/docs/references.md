# Algorithm References

FastNSA implements and adapts algorithms described in the publications below.
The papers define the algorithmic ideas and complexity results; the FastNSA API,
validation, result normalization, cancellation, and persistence behavior are
documented by this package and may differ from research artifact interfaces.

When reporting scientific results, cite the paper corresponding to the selected
method and record the FastNSA version, parameters, and random seed.

## Point-to-Line Nearest-Edge Matching

The point-to-line nearest-neighbor matching methods are based on:

Hongwei Ye, Tsz Nam Chan, Dingming Wu, Leong Hou U, and Ruisheng Wang. 2026.
**An Efficient and Accurate Grid Compression Solution for Point-to-Line Nearest
Neighbor Search.** In *Proceedings of KDD 2026*, 6152-6163.
[DOI: 10.1145/3770855.3817711](https://doi.org/10.1145/3770855.3817711).

Research artifact:
[fastnsa/grid-comp-plnns](https://github.com/fastnsa/grid-comp-plnns).

Relevant API: MappingMethod, NativePointSet.map_xy, PointEvents, and
preprocess_events.

## Network Kernel Density Visualization

ADA, IA, and their combined augmentation family are described in:

Tsz Nam Chan, Zhe Li, Leong Hou U, Jianliang Xu, and Reynold Cheng. 2021.
**Fast Augmentation Algorithms for Network Kernel Density Visualization.**
*PVLDB* 14(9), 1503-1516.
[DOI: 10.14778/3461535.3461540](https://doi.org/10.14778/3461535.3461540).

LION is described in:

Tsz Nam Chan, Rui Zang, Bojian Zhu, Leong Hou U, Dingming Wu, and Jianliang Xu.
2024. **LION: Fast and High-Resolution Network Kernel Density Visualization.**
*PVLDB* 17(6), 1255-1268.
[DOI: 10.14778/3648160.3648168](https://doi.org/10.14778/3648160.3648168).

Research artifact: [edisonchan2013928/LION](https://github.com/edisonchan2013928/LION).

Relevant API: NKDVMethod and network_kde. FastNSA currently exposes only the
validated Epanechnikov kernel across optimized methods.

## Network K-function

The NS optimization, ASPS framework, and network K-function baselines are
described in:

Tsz Nam Chan, Leong Hou U, Yun Peng, Byron Choi, and Jianliang Xu. 2022.
**Fast Network K-function-based Spatial Analysis.** *PVLDB* 15(11), 2853-2866.
[DOI: 10.14778/3551793.3551836](https://doi.org/10.14778/3551793.3551836).

Research artifact:
[edisonchan2013928/Network-K-function](https://github.com/edisonchan2013928/Network-K-function).

The sampling and threshold-sharing method is described in:

Hongwei Ye, Tsz Nam Chan, Leong Hou U, Dingming Wu, Wei Tu, Ruisheng Wang, and
Joshua Zhexue Huang. 2026. **STARS: A Sampling and Threshold Sharing Solution
for Network K-function Analytics.** In *Proceedings of KDD 2026*.
[DOI: 10.1145/3770855.3817709](https://doi.org/10.1145/3770855.3817709).

Research artifact: [fastnsa/STARS](https://github.com/fastnsa/STARS).

Relevant API: KFunctionMethod and network_k_function. STARS is approximate;
report epsilon together with the method and FastNSA version.

## Spatiotemporal Network K-function

EAR, MTS, and their combined EARTH method are described in:

Tsz Nam Chan, Hongwei Ye, Leong Hou U, Yun Peng, Dingming Wu, Jianliang Xu, and
Christian S. Jensen. 2026. **EARTH: Accelerating Spatiotemporal Network
K-function-based Analytics.** In *Proceedings of KDD 2026*, 297-308.
[DOI: 10.1145/3770855.3817706](https://doi.org/10.1145/3770855.3817706).

Research artifact: [edisonchan2013928/EARTH](https://github.com/edisonchan2013928/EARTH).

Relevant API: SpatiotemporalKFunctionMethod and
spatiotemporal_network_k_function. RQS and SPS are baseline methods discussed
by the paper; MTS and EARTH require regular spatial and temporal threshold axes
in the FastNSA implementation.

## Citing FastNSA

Algorithm papers should be cited for their methods, but they do not replace a
software citation. Until a dedicated FastNSA archival citation is published,
also report the package name, version, repository URL, and access date. A future
release may provide CITATION.cff and a version-specific archival DOI.

