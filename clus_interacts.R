#clus_interacts() uses SingleCellSignalIR to identify ligand-receptor
#pairs among sc clusters
clus_interact <- function(matrix, genelist, clusterlist, ...) {
  # Ligand/Receptor analysis using SingleCellSignalR
  signal = cell_signaling(data=matrix,genes=genelist,cluster=clusterlist)
  inter.net <- inter_network(data = matrix, signal = signal, genes = genelist, cluster = clusterlist, write = FALSE)
  # Visualization
  visualize_interactions(signal = signal)
  datalist <- list(signal, inter.net)   # return value 
  return(datalist)
}
