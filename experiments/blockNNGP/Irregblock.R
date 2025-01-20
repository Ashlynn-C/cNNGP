
#########################################################
#***** functions for irregular blocks 
#***** with similar number of points on each block
#***** treenew has the kd-tree of the set of locations
#***** for now you can only set 
#***** n.blocks = 2^2, 2^3,2^4,2^5,2^6,2^7 
#########################################################

part_v = function(treenew,level_k,num_block,blocks){
  
  k                 <- level_k[1]
  vb                <- treenew[k,]$x
  indblock          <- which(loc[,1] <= vb)
  blocks[indblock]  <- num_block [1]
  blocks[-indblock] <- num_block [4]
  
  k                 <- level_k[2]
  hb                <- treenew[k,]$y
  indblock1         <- which(loc[,2] > hb)
  indblock2         <- which(blocks==num_block [1])
  indblock          <- intersect(indblock1, indblock2)
  blocks[indblock]  <- num_block [2]
  
  k                 <- level_k[3]
  hb                <- treenew[k,]$y
  indblock1         <- which(loc[,2] < hb)
  indblock2         <- which(blocks==num_block [4])
  indblock          <- intersect(indblock1, indblock2)
  blocks[indblock]  <- num_block [3]
  return(blocks)
  
}

part_v2 = function(treenew,level_k,num_block,blocks){
  
  k                 <- level_k[1]
  vb		  <- treenew[k,]$x
  indblock1 	  <- which(loc[,1] > vb)
  indblock2 	  <- which(blocks==num_block[1])
  indblock 	  <- intersect(indblock1, indblock2)
  blocks[indblock]  <- num_block[2]
  
  return(blocks)
  
}


part_h2 = function(treenew,level_k,num_block,blocks){
  
  k 		 <- level_k[1]
  hb	         <- treenew[k,]$y
  indblock1 	 <- which(loc[,2] > hb)
  indblock2 	 <- which(blocks==num_block[1])
  indblock 	 <- intersect(indblock1, indblock2)
  blocks[indblock] <- num_block[2]
  return(blocks)
  
}


kdtree_blocks =function(treenew,n.blocks,loc){
  
  nexp = log(n.blocks)/log(2)
  
  # 4 blocks
  m = n.blocks/4
  blocks 		 <- matrix(NA,n,1)
  level_k	 	 <- c(1,2,3)
  num_block 	 <- c(1,m+1,2*m+1,3*m+1)
  blocks 		 <- part_v(treenew,level_k,num_block,blocks)
  sort(unique(blocks))
  
  if(nexp >= 3){
    # 8 blocks
    m		<-4
    k 		<- n.blocks/4
    r 		<- k-1
    if(n.blocks == 8)   t=1
    if(n.blocks == 16)  t=1
    if(n.blocks == 32)  t=4
    if(n.blocks == 64)  t=8
    if(n.blocks == 128) t=16
    for (i in 1:m){
      level_k 	<- 3+i 
      num_block 	<- c((k*i-r),(k*i-r+t))
      blocks		<- part_v2(treenew,level_k,num_block,blocks)
    }
    sort(unique(blocks))
  }
  
  if(nexp >= 4){
    # 16 blocks
    if(n.blocks == 16)  t=2
    if(n.blocks == 32)  t=2
    if(n.blocks == 64)  t=4
    if(n.blocks == 128) t=8
    blocks1 	<- sort(unique(blocks))
    m 		<- 8
    for (i in 1:m){
      level_k 	<-  c(7+i)
      num_block 	<- c(blocks1[i],blocks1[i]+t)
      blocks		<- part_h2(treenew,level_k,num_block,blocks)
    }
    sort(unique(blocks))
  }
  
  
  if(nexp >= 5){
    # 32 blocks 
    if(n.blocks == 32)  t=1
    if(n.blocks == 64)  t=1
    if(n.blocks == 128) t=4
    blocks1 	<- sort(unique(blocks))
    m 		<- 16
    for (i in 1:m){
      level_k 	<- c(15+i)
      num_block 	<- c(blocks1[i],blocks1[i]+t)
      blocks		<- part_v2(treenew,level_k,num_block,blocks)
    }
    sort(unique(blocks))
  }
  
  
  if(nexp >= 6){
    # 64 blocks 
    if(n.blocks == 64)  t=2
    if(n.blocks == 128) t=2
    blocks1 	<- sort(unique(blocks))
    m 		<- 32
    for (i in 1:m){
      level_k 	<- c(31+i)
      num_block 	<- c(blocks1[i],blocks1[i]+t)
      blocks		<- part_h2(treenew,level_k,num_block,blocks)
    }
    sort(unique(blocks))
  }
  
  
  if(nexp >= 7){
    # 128 blocks 
    if(n.blocks == 128) t=1
    blocks1 	<- sort(unique(blocks))
    m 		<- 64
    for (i in 1:m){
      level_k 	<- c(63+i)
      num_block 	<- c(blocks1[i],blocks1[i]+t)
      blocks		<- part_v2(treenew,level_k,num_block,blocks)
    }
    sort(unique(blocks))
  }
  
  return(blocks)
}

#' k-d tree
#' function from package Mathart in R
#'
#' Computes a k-d tree for a given set of points
#' @param points A data frame with columns for x and y coordinates, and each point in a row. Refer to the \href{https://en.wikipedia.org/wiki/K-d_tree}{Wikipedia article} for details.
#' @keywords k-d tree
#' @export
#' @examples
#' kdtree()

kdtree = function(points, minmax = FALSE) {
  n 	<- nrow(points)
  df 	<- data.frame(xmin = numeric(n), xmax = numeric(n), 
                    ymin = numeric(n), ymax = numeric(n),
                    dir  = character(n),
                    x    = numeric(n), y    = numeric(n), 
                    xend = numeric(n), yend = numeric(n),
                    minmax = character(n)) %>%
    mutate(dir = as.character(dir), minmax = as.character(n))
  
  i <- 1
  k <- 2
  l <- n
  
  df[1, c("xmin", "xmax", "ymin", "ymax")] <- c(min(points$x), max(points$x),
                                                min(points$y), max(points$y))
  df[1, "dir"] <- "v"
  df[1, "minmax"] <- "min"
  
  while(i < l) {
    if(df$dir[i] == "v") {
      temp <- points %>%
        filter(x > df[i, "xmin"], x < df[i, "xmax"], y > df[i, "ymin"], y < df[i, "ymax"])
      if(minmax & df$minmax[i] == "min") {
        temp <- temp %>% summarise(x = min(x) + 1)
      } else if (minmax & df$minmax[i] == "max") {
        temp <- temp %>% summarise(x = max(x) - 1)
      } else {
        temp <- temp %>% summarise(x = median(x))
      }
      df[i, c("x", "xend", "y", "yend")] <- c(temp, temp, df$ymin[i], df$ymax[i])
      df[k, ] 		<- df[i, ]
      df[k+1, ] 	<- df[i, ]
      df$xmax[k] 	<- temp
      df$xmin[k+1] 	<- temp
      df$dir[k] 	<- "h"
      df$dir[k+1] 	<- "h"
      df$minmax[k] 	<- ifelse(df$minmax[i] == "min", "max", "min")
      df$minmax[k+1] 	<- ifelse(df$minmax[i] == "min", "min", "max")
      k <- k + 2
      i <- i + 1
    } else {
      temp <- points %>%
        filter(x > df[i, "xmin"], x < df[i, "xmax"], y > df[i, "ymin"], y < df[i, "ymax"])
      if(minmax & df$minmax[i] == "min") {
        temp <- temp %>% summarise(y = min(y) + 1)
      } else if (minmax & df$minmax[i] == "max") {
        temp <- temp %>% summarise(y = max(y) - 1)
      } else {
        temp <- temp %>% summarise(y = median(y))
      }
      df[i, c("x", "xend", "y", "yend")] <- c(df$xmin[i], df$xmax[i], temp, temp)
      df[k, ] 	     <- df[i, ]
      df[k+1, ]      <- df[i, ]
      df$ymax[k]     <- temp
      df$ymin[k+1]   <- temp
      df$dir[k]      <- "v"
      df$dir[k+1]    <- "v"
      df$minmax[k]   <- ifelse(df$minmax[i] == "min", "max", "min")
      df$minmax[k+1] <- ifelse(df$minmax[i] == "min", "min", "max")
      k <- k + 2
      i <- i + 1
    }
  }
  df
}
###################################################


