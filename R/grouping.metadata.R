# Functions to deduce metadata on groupings in dataframe.

grouping.metadata_nonest <- function(boostdf){
  groups.metadata <- unique(boostdf[,"id"])
  I.metadata <- length(groups.metadata)
  n_i.metadata <- rep(0,I.metadata)
  for (i in seq_len(I.metadata)) {
    n_i.metadata[i] <- nrow(boostdf[boostdf$id==groups.metadata[i],])
  }
  n_obs.metadata <- sum(n_i.metadata)
  return(list(groups.metadata=groups.metadata, I.metadata=I.metadata, n_i.metadata=n_i.metadata, n_obs.metadata=n_obs.metadata))
}

grouping.metadata_nested <- function(boostdf){
  groups.metadata <- unique(boostdf[,"id"])
  I.metadata <- length(groups.metadata)
  n_ij_limits <- length(unique(boostdf[,"subid"]))
  J.metadata <- rep(0,I.metadata)
  n_ij.metadata <- matrix(nrow=I.metadata,ncol=n_ij_limits)
  for (i in seq_len(I.metadata)) {
    subboostdf <- boostdf[boostdf$id==groups.metadata[i],]
    subgroups.metadata <- unique(subboostdf[,"subid"])
    J.metadata[i] <- length(subgroups.metadata)
    for (j in seq_len(J.metadata[i])) {
      n_ij.metadata[i,j] <- nrow(subboostdf[subboostdf$subid==subgroups.metadata[j],])
    }
  }
  n_obs.metadata <- sum(n_ij.metadata, na.rm=TRUE)
  return(list(groups.metadata=groups.metadata, I.metadata=J.metadata, n_i.metadata=n_ij.metadata, n_obs.metadata=n_obs.metadata))
}
