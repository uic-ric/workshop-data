overlap_index <- function( cluster1, cluster2 ){
   # make factor vectors from clusters
   c1 = as.factor(cluster1)
   c2 = as.factor(cluster2)
   # define a set of "names" for our cells, which we'll compare between clusters
   # the names are arbitraty, so we'll just number them
   names = 1:length(c1)
   # make distance matrix to store our calculations in
   my_dist = matrix(nrow=length(levels(c1)),ncol=length(levels(c2)))
   rownames(my_dist) = levels(c1)
   colnames(my_dist) = levels(c2)
   for(i in 1:length(levels(c1))){
      for(j in 1:length(levels(c2))){
         # get the list of cells in each cluster
         c1.sub = names[c1==levels(c1)[i]]
         c2.sub = names[c2==levels(c2)[j]]
         # get the intersection and min cluster size
         int = length(intersect(c1.sub,c2.sub))
         minsize = min(length(c1.sub),length(c2.sub))
         # overlap index
         overlap = int/minsize
         # store in matrix
         my_dist[i,j] = overlap
      }
   }
   return(my_dist)
}

