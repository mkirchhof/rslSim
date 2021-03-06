# pRSL: Approximate Inference Simulation Experiments

This repo contains the code to run experiments to benchmark the approximation quality of the noisy-or loopy belief propagation algorithm used in pRSL. Each folder is a self-contained experiment, that is it includes a (possibly old) version of rsl.R and norn.R, an experiment.R file that runs the actual experiment and generates simRes.RData files as results and an analysis.R file that summarizes those results. Additionally, there is the generateDatasetsNoisyOR.R which generates the datasets stored in the /small, /medium and /almostlarge folder, as detailed in the paper.

So, if you want to reproduce the experiments from scratch:
1. run generateDatasetsNoisyOR.R
2. Choose an experiment folder
3. Run experiment.R
4. Run analysis.R

# Further Resources:

- Main pRSL code repo (with up-to-date version of pRSL): https://github.com/mkirchhof/rsl

- Approximate inference experiments: https://github.com/mkirchhof/rslSim

- Benchmarks of pRSL on multi-label datasets: https://github.com/mkirchhof/rslBench

- Application of pRSL for transfer learning in human activity recognition: https://github.com/mkirchhof/rslAppl

- Paper: Kirchhof, M., Schmid, L., Reining, C., ten Hompel, M., Pauly, M.: pRSL: Interpretable Multi–label Stacking by Learning Probabilistic Rules. In: Proceedings of the 37th Conference on Uncertainty in Artificial Intelligence, PMLR, in press (2021).
