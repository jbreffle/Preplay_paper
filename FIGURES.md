# Figures

A description of how each figure in the manuscript is generated.

## Figure 1: Illustration of the randomly clustered model

- a,b,d-g: `plotNetSchematics.m`
- c: Created directly in Inkscape
- h: `swiGrid.m` with option `panelToPlot='1g'`

## Figure 2: Spatially correlated reactivations in networks without environment specific connectivity or plasticity

- a,b: `simulateNetwork.m` with `plotFig2=true`
  - Vm from cells 131 and 275
- c-f: `gridAnalysisCombNet.m` with simulation grid data "mainPaperGrid"
  - `paramSetInds = combvec([2], [3])'`;
  - d,e: Event from network 4 (set `plotNetsBest = true`)

## Figure 3: The model produces place fields with similar properties to hippocampal place fields

- a: `exptAnalysisPlaceFields.m` with epoch 2
- c-b: `gridAnalysisPlaceFieldStats.m` with simulation grid data "mainGrid300s"
  - b is parameter point `analysisParam.paramPoint = [2, 3]`

## Figure 4: Preplay depends on modest cluster overlap

- a-b: `exptSleepDecodePlotting.m` with `analysisParam.trajToInclude = 1:4`
- c-d: `gridAnalysisCombNet.m` with simulation grid data "mainGrid300s"
  - `paramSetInds = combvec([2], [3])';`
- e-f: `gridAnalysisDecoding.m` with simulation grid data "mainGrid300s"

## Figure 5: Coherent spiking within clusters supports preplay

- `gridAnalysisEventSpiking.m` simulation grid data "mainPaperGrid"
- a: First event from parameter point [2,3] with `plotAndPauseOnExample = true`
- b,c: From parameter point [2,3]
- c,e: From all parameter points

## Figure 6: Preplay is abolished when events are decoded with shuffled cell identities but is preserved if cell identities are shuffled only within clusters

- All plotted by `plotShuffledDecoding.m`
- a: `shuffleToPlot = "clusterIndependent"`
- b: `shuffleToPlot = "withinCluster"`
- c: `shuffleToPlot = "singleClusterCells"`

## Figure 7: Place cells’ mean event rank are correlated with their place field location when accounting for decode direction

- a-b: Can be plotted with `gridAnalysisEventSpikingSequences.m` with simulation grid data "mainGrid300s"
  - a: analysisParam.invertDirectionality=false
  - b: analysisParam.invertDirectionality=true
- c-d: `grid_2023-07-22T14-50shuffledDecode2024-04-26T21-15` plotted by `plotShuffledDecodingSpikeSequences.m`

## Figure 8: The Small-World Index of networks correlates with preplay quality

- a,b,c left column: `swiGrid.m`
- a: `gridAnalysisDecoding.m` with simulation grid data "mainGridConnProb04"
- b: `gridAnalysisDecoding.m` with simulation grid data "mainPaperGrid"
- c: `gridAnalysisDecoding.m` with simulation grid data "mainGridConnProb12"

## Figure 9: Trajectories decoded from population-burst events are significantly correlated with linear trajectories in arbitrary environments

- All from simulation grid data "mainPaperGrid"
- a,b,c,e: `gridAnalysisCombNet.m` with `paramSetInds = combvec([2], [3])'`
  - a,b: Use `plotNetPFs = true`
  - a: Network 3
  - c: `plotNetsBest = true`, network 2, delete title in Inkscape
  - e: `figSettings = 'manuscriptAlt'`, and run with each of `ithEnv` equal to 1, 2, 3, and 4
- d: `gridAnalysisDecoding.m` with `analysisParam.paramPoint = [2, 3]`

## Supplemental figures

## Figure 1--figure supplement 1: Comparison of the randomly clustered network and the canonical Watts-Strogatz small-world network

- All panels created by running `scripts/misc-code/plotMiniNetworks.m`

## Figure 3-—figure supplement 1: The simulated cells have greater place information than time information

- All panels created by running `scripts/variedPFSpeed.m`

## Figure 4—figure supplement 1: Example preplay events from the Shin et al., 2019 data

- Selected events from running `scripts/sleep-decoding/exptSleepDecodePlotting.m`

## Figure 4--figure supplement 2: Significant preplay can typically be identified with as 1082 few as 50 cells

- The script `scripts/varyDecodeParam.m` generates the data for this figure
- The script `scripts/plotVariedDecodeParam.m` plots the figure

## Figure 4—figure supplement 3: Preplay statistics by trajectory for Shin et al., 2019 data

- Each row of panels is plotted by running `scripts/sleep-decoding/exptSleepDecodePlotting.m` with `analysisParam.trajToInclude` set to 1, 2, 3, or 4 for each row of panels

## Figure 4--figure supplement 4: Additional simulations support the consistency and robustness of the model to variations in spatial input forms

- Input schematics produced by `plotNetSchematics.m`
- Parameter grid produced by `gridAnalysisDecoding.m`
- CDF and p-val matrix produced by `gridAnalysisCombNet.m`
- a: `mainGridLinearTernary`
- b: `mainGridStepped5`
- c: `mainGridStepped3`
- d: `mainGridStepped1`
- e: `mainGridNoClusterBias`

## Figure 5--figure supplement 1: Relationship between cluster activation and preplay

- All panels created by running `gridAnalysisEventSpiking.m` with simulation grid data "mainGrid300s"
