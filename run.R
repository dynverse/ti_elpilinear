#!/usr/local/bin/Rscript

task <- dyncli::main()

library(ElPiGraph.R)
library(dplyr)
library(purrr)
library(readr)


#   ____________________________________________________________________________
#   Load data                                                               ####

expression <- as.matrix(task$expression)
parameters <- task$parameters

checkpoints <- list()
checkpoints$method_afterpreproc <- as.numeric(Sys.time())

#   ____________________________________________________________________________
#   Infer the trajectory                                                    ####
principal_graph_function <- computeElasticPrincipalCurve

# infer the principal graph, from https://github.com/Albluca/ElPiGraph.R/blob/master/guides/base.md
principal_graph <- principal_graph_function(
  X = expression,
  NumNodes = parameters$NumNodes,
  NumEdges = parameters$NumEdges,
  InitNodes = parameters$InitNodes,
  MaxNumberOfIterations = parameters$MaxNumberOfIterations,
  eps = parameters$eps,
  CenterData = parameters$CenterData,
  Lambda = parameters$Lambda,
  Mu = parameters$Mu,
  drawAccuracyComplexity = FALSE,
  drawEnergy = FALSE,
  drawPCAView = FALSE,
  n.cores = 1
)

# compute pseudotime, from https://github.com/Albluca/ElPiGraph.R/blob/master/guides/pseudo.md
PartStruct <- PartitionData(
  X = expression,
  NodePositions = principal_graph[[1]]$NodePositions,
  nCores = 1
)

ProjStruct <- project_point_onto_graph(
  X = expression,
  NodePositions = principal_graph[[1]]$NodePositions,
  Edges = principal_graph[[1]]$Edges$Edges,
  Partition = PartStruct$Partition
)

checkpoints$method_aftermethod <- as.numeric(Sys.time())

#   ____________________________________________________________________________
#   Process & save the model                                               ####
milestone_network <- ProjStruct$Edges %>%
  as_tibble() %>%
  rename(from = row, to = col) %>%
  mutate(
    from = paste0("M", from),
    to = paste0("M", to),
    length = ProjStruct$EdgeLen,
    directed = FALSE
  )

progressions <- tibble(cell_id = rownames(expression), edge_id = ProjStruct$EdgeID) %>%
  left_join(milestone_network %>% select(from, to) %>% mutate(edge_id = row_number()), "edge_id") %>%
  select(-edge_id) %>%
  mutate(percentage = pmin(1, pmax(0, ProjStruct$ProjectionValues)))

output <- lst(
  cell_ids = rownames(expression),
  milestone_network,
  progressions,
  timings = checkpoints
)

dynwrap::wrap_data(cell_ids = rownames(expression)) %>%
  dynwrap::add_trajectory(
    milestone_network = milestone_network,
    progressions = progressions
  ) %>%
  dynwrap::add_timings(timings = checkpoints) %>%
  dyncli::write_output(task$output)
