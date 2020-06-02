# For more details on graph visualization see tutorial here:
# https://kateto.net/network-visualization

library('igraph')

pairs_df <- data.frame((pairs))
net <- graph_from_data_frame(d=pairs_df, directed=T) 
net <- simplify(net, remove.multiple = F, remove.loops = T)

node_counts <-data.frame(table(pairs_df$X1))
pal1 <- heat.colors(max(node_counts$Freq)), alpha=1)
V(net)$color <- pal1[node_counts$Freq]

l <- layout_with_graphopt(net, charge=0.1)
plot(net, edge.arrow.size=.2, edge.width=1, edge.curved=.1,vertex.label=NA, 
     vertex.size=6, margin=0, layout=layout_with_mds, col=pal1)

deg.dist <- degree_distribution(net, cumulative=T, mode="all")
plot( x=0:max(degree(net)), y=1-deg.dist, pch=19, cex=1.2, col="orange", 
      xlab="Degree", ylab="Cumulative Frequency")


library('visNetwork') 
visNetwork(pairs_df, pairs_df, width="100%", height="400px")
