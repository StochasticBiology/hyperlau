#' Reconstruct an ancestral state based on two specified daughters
#' 
#' Data is specified as binary strings of equal length, which may include "?" for uncertainty
#'
#' @param s1 A binary string, possibly with "?"s
#' @param s2 A binary string, possibly with "?"s
#'
#' @return A binary string describing the ancestor state
#' @examples
#' make.ancestor("101", "0??")
#' @export
make.ancestor = function(s1, s2) {
  ss1 = strsplit(s1, split="")[[1]]
  ss2 = strsplit(s2, split="")[[1]]
  news = rep("0", length(ss1))
  # truth table
  # s1  000111???
  # s2  01?01?01?
  # new 00001?0??
  for(i in 1:length(ss1)) {
    if(ss1[i]=="1") {
      if(ss2[i]=="1") { news[i] = "1" }
      if(ss2[i]=="?") { news[i] = "?" }
    } else if(ss2[i] == "1" && ss1[i] == "?") { 
      news[i] = "?" 
    } else if(ss1[i]=="?" && ss2[i] == "?") { news[i] = "?" }
  }
  return(paste(news, collapse=""))
}

#' Extract transitions from a dataset given a phylogeny and (possibly uncertain) observations on its tips
#' 
#' @param data A dataframe. The first column should give a label that matches the tip labels of the phylogeny. The second should be a character string containing "0"s, "1"s, and "?"s
#' @param tree A tree linking the observations, with tip labels corresponding to the first column of the data frame
#'
#' @return A named list containing independent transitions, the source data, and reconstructed node states 
#' @export
curate.uncertain.tree = function(data, tree) {
  my.data = data
  my.rooted.tree = tree
  colnames(my.data)[1] = "label"
  # check for missing data in both barcodes and tips... from curate.tree
  match.set = match(my.data$label, my.rooted.tree$tip.label)
  if(any(is.na(match.set))) {
    message("Found observations that didn't correspond to tips of the tree!")
    my.data = my.data[-which(is.na(match.set)),]
    match.set = match(my.data$label, my.rooted.tree$tip.label)
  }
  
  if(any(duplicated(my.data$label))) {
    message("Duplicates in observation set!")
  }
  
  # prune tree to include only those tips in the barcode dataset
  my.tree = ape::drop.tip(my.rooted.tree,
                          my.rooted.tree$tip.label[-match.set])
  
  final.df = my.data
  
  # assign arbitrary labels to nodes for bookkeeping
  my.tree$node.label = as.character(1:my.tree$Nnode)
  tree.labels = c(my.tree$tip.label, my.tree$node.label)
  tree.states = rep("", length(tree.labels))
  
  ##########
  ## assign states (with uncertainty) to tips and then nodes of tree
  
  # populate tree states with states known from tip data
  df = data.frame()
  for(i in 1:nrow(final.df)) {
    this.line = final.df[i,]
    df = rbind(df, data.frame(label=this.line[1], state=this.line[length(this.line)]))
    ref = which(tree.labels == as.character(this.line[1]))
    if(length(ref) > 0) {
      tree.states[ref] = as.character(this.line[2])
    }  
  }
  
  # first, reconstruct full tree
  trans.df = data.frame()
  # assume we've got work to do reconstructing ancestors
  change = TRUE
  # loop until we're done
  while(change == TRUE) {
    change = FALSE
    # which nodes don't yet have states
    todo = which(tree.states == "")
    if(length(todo) == 0) {
      change = FALSE
      break
    }
    # loop through empty nodes
    for(i in 1:length(todo)) {
      # get children of this node, and figure out if their states are specified
      ref = todo[i]
      kids = phangorn::Children(my.tree, ref)
      kid.empty = which(tree.states[kids] == "")
      if(length(kids) == 2 && length(kid.empty) == 0) {
        # if the childrens' states are specified, reconstruct this ancestral state and add transitions to the dataframe
        print(paste("Tree", 1, ": populating", ref, "based on kids", kids[1], "and", kids[2], collapse = " "))
        tree.states[ref] = make.ancestor(tree.states[kids[1]], tree.states[kids[2]])
        trans.df = rbind(trans.df, data.frame(tree=1, b = ref, a = kids[1], before=tree.states[ref], after=tree.states[kids[1]]))
        trans.df = rbind(trans.df, data.frame(tree=1, b = ref, a = kids[2], before=tree.states[ref], after=tree.states[kids[2]]))
        change = TRUE
      }
    }
  }
  
  ###########
  ### pull set of independent transitions 
  
  # now, consider which transitions we need to drop from this set to make sure we're dealing with indepedent behaviours
  
  # first, consider only those transitions to a tip state from its immediate ancestor...
  prune.df = trans.df[trans.df$a < length(my.tree$tip.label),]
  # ... and those that actually contain information (i.e. before and after states are different)
  prune.df = prune.df[prune.df$before != prune.df$after,]
  
  unsafe = prune.df
  pairs = data.frame()
  # go through every remaining transition
  for(i in 1:nrow(unsafe)) {
    this.safe = TRUE
    if(nrow(pairs) > 0) {
      # if the ancestral *node* (not state) is already in our working set, ignore this transition
      # (avoid double counting P(0->a))
      if(unsafe$b[i] %in% pairs$b) {
        this.safe = FALSE
      }
      # go through transitions in our working set, comparing current transition with each
      for(j in 1:nrow(pairs)) {
        a1 = unsafe$b[i]
        a2 = pairs$b[j]

        # if the state of the ancestors' ancestor isn't compatible with 0^L, ignore this transition
        a12 = ape::getMRCA(my.tree, c(a1, a2))
        
        if(grepl("1|\\?", tree.states[a12])) {
          this.safe = FALSE
        }
      }
    }
    # if we're not coupled to any existing transitions in the working set, add to working set
    if(this.safe == TRUE) {
      pairs = rbind(pairs, unsafe[i,])
    }
  }
  
  # confirm that we're only looking at ancestors of transition pairs that are compatible with 0^L
  a.s = c()
  ref = 2
  for(i in 1:nrow(pairs)) {
    if(i != ref) {
      a.s = c(a.s, tree.states[ape::getMRCA(my.tree, c(pairs$b[ref], pairs$b[i]))])
    }
  }
  
  list.r = list(independents = pairs[,c(4,5)],
                tree = my.rooted.tree,
                data = my.data,
                states= tree.states)
  return(list.r)
}

