# draw tree

drawTree <- function(tree, xlim=NULL, mark=NULL) {
  
  g1 <- as(tree, 'phylo4')
  d  <- data.frame(color=rep("black", length(tree$tip.label)), row.names=tree$tip.label)
  d$color <- as.character(d$color)
  
  if(!is.null(mark)) d$color[grepl(mark, rownames(d))] <- "red" # annotate tip
  
  g2 <- phylo4d(g1, d)
  rNodeData <- data.frame(color=rep("black", nNodes(g1)), row.names=nodeId(g1, "internal"))
  rNodeData$color <- as.character(rNodeData$color)
  
  if(!is.null(mark)) rNodeData$color[which(rownames(rNodeData) %in% names(nodeLabels(g1))[grepl(mark, nodeLabels(g1))])] <- "red" # annotate node
  nodeData(g2) <- rNodeData
  
  print(g2)
  
  p <- ggtree(g2, aes(color=I(color)), size=0.3) #+ geom_nodelab() 
  
  if(!is.null(xlim)) p <- p + scale_x_continuous(limits=c(0,xlim), oob=scales::squish)
  
  edge           <- data.frame(tree$edge, edge_length=tree$edge.length)
  colnames(edge) <- c("parent", "node", "edge_length")
  #p              <- p %<+% edge + geom_label(aes(x=branch, label=edge_length), na.rm = TRUE)
  
  #print(p)
  
  return(p)
  
}
