#count_signal will count the number of interactions for a given cluster with ever other cluster in a sinle cell data set
count_signal <- function(signallist, clus_no) {
  signal_count <- 0
  for (cell_inter in signallist) {
    if (length(grep(paste0("cluster ",clus_no,"-"), cell_inter))>0) {
      signal_count <- signal_count +1
    }
  }
  return(signal_count)
}
