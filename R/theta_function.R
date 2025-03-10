library(ape)
library(phytools)
library(expm)

#'Theta and zeta indices for phylogenetic signal
#'
#'@param cat A dataframe containing data for the categorical trait. The data preferably should be in different columns for each class of the trait and 0 indicating the species does not have this class, while 1 that it does. A species could have more than one categories.
#'@param tree The phylogenetic tree in phylo format. The tree must include all species in the trait dataset
#'@param matching The similarity criterion to be used "complete" for assuming two species as identical if they match in all classes, and "partial" if we assume that species are considered similar if they share at least one class of the trait
#'
#'@return The values of the theta and zeta indices. The values are between 0 and 1, and the smaller the value the stronger the signal
#'
#'@examples
#'tree <- rtree(50)
#'cat <- matrix(sample(c(0,1), 200, replace=TRUE), ncol=4, nrow=50)
#'cat <-data.frame(cat)
#'rownames(cat)<-tree$tip.label
#'theta(cat, tree, "complete")
#'theta(cat, tree, "partial")

theta<-function(cat, tree, matching="complete"){

if(matching!="partial" & matching!="complete"){
stop("Matching criterion for trait comparison among species could be only partial or complete")
}

if(is.data.frame(cat)==FALSE){
stop("Trait data need to be in data frame format, with row names being the species names")
}

#remove lines from trait dataset that have NA
cat<-na.omit(cat)

if(sum(row.names(cat) %in% tree$tip.label)/dim(cat)[1] < 1){
stop("The list of species in trait dataset is not present in tree")
}

if(sum(row.names(cat) %in% tree$tip.label)/length(tree$tip.label) < 1){
tree<-drop.tip(tree,tree$tip.label[-match(row.names(cat), tree$tip.label)])
}


original_tree <- tree


##### merge tips with same value for trait, iteratively to minimum tree

repeat{   # the do while loop starts here

is_terminal <- sapply(tree$edge[,2], function(x) !(x %in% tree$edge[,1]))
terminal_nodes <- tree$edge[is_terminal,1]
count<-1
terminals<-numeric()
lg<-length(terminal_nodes)
t2<-terminal_nodes[2:lg]==terminal_nodes[1:(lg-1)]
count=sum(t2)
terminals<-terminal_nodes[t2]
terminals<-terminals[1:count]

lg<-count
count<-1
same<-numeric()

for(j in 1:lg){
test<-tree$edge[,1]==terminals[j]

if(matching=="partial"){
 if(sum(cat[tree$tip.label[tree$edge[test,2][1]],]*cat[tree$tip.label[tree$edge[test,2][2]],])>0){
   same[count] <- terminals[j]
   count <- count +1
   newvalue <- cat[tree$tip.label[tree$edge[test,2][1]],]*cat[tree$tip.label[tree$edge[test,2][2]],]
   cat[tree$tip.label[tree$edge[test,2][1]],] <- newvalue
   cat[tree$tip.label[tree$edge[test,2][2]],] <- newvalue
 }
}
else{
if(matching=="complete"){
 if(sum(cat[tree$tip.label[tree$edge[test,2][1]],]==cat[tree$tip.label[tree$edge[test,2][2]],])==length(cat[tree$tip.label[tree$edge[test,2][2]],])){
   same[count] <- terminals[j]
   count <- count +1
 }
}
}  #end else of comparisons

}


new_tree<-tree
droping<-numeric()
if(length(same)>0){
 for(k in 1:length(same)){
   test<-tree$edge[,1]==same[k]
   a<-max(test*tree$edge.length)
   tree$edge.length[test]<-a
   droping[k]<-tree$edge[test,2][1]
 }
   new_tree<-drop.tip(tree, droping)
}


if(length(new_tree$tip.label) == length(tree$tip.label)){break} #ending the do while loop
else{tree <-new_tree}
}   # end of do while


theta<-sum(tree$edge.length) / sum(original_tree$edge.length)
zeta<-(dim(tree$edge)[1]+1)/(dim(original_tree$edge)[1]+1)

output <- data.frame(theta, zeta)
return(output)
}



###########################################################################
###########################################################################

## lighter version for estimating theta and zeta for randomizations
## without data checking
theta2<-function(species, trait, tree, matching){

original_tree <- tree
cat<-data.frame(trait)
rownames(cat)<-species

repeat{   # the do while loop starts here

is_terminal <- sapply(tree$edge[,2], function(x) !(x %in% tree$edge[,1]))
terminal_nodes <- tree$edge[is_terminal,1]
count<-1
terminals<-numeric()
lg<-length(terminal_nodes)
t2<-terminal_nodes[2:lg]==terminal_nodes[1:(lg-1)]
count=sum(t2)
terminals<-terminal_nodes[t2]
terminals<-terminals[1:count]

lg<-count
count<-1
same<-numeric()

for(j in 1:lg){
test<-tree$edge[,1]==terminals[j]

if(matching=="partial"){
 if(sum(cat[tree$tip.label[tree$edge[test,2][1]],]*cat[tree$tip.label[tree$edge[test,2][2]],])>0){
   same[count] <- terminals[j]
   count <- count +1
   newvalue <- cat[tree$tip.label[tree$edge[test,2][1]],]*cat[tree$tip.label[tree$edge[test,2][2]],]
   cat[tree$tip.label[tree$edge[test,2][1]],] <- newvalue
   cat[tree$tip.label[tree$edge[test,2][2]],] <- newvalue
 }
}
else{
if(matching=="complete"){
 if(sum(cat[tree$tip.label[tree$edge[test,2][1]],]==cat[tree$tip.label[tree$edge[test,2][2]],])==length(cat[tree$tip.label[tree$edge[test,2][2]],])){
   same[count] <- terminals[j]
   count <- count +1
 }
}
}  #end else of comparisons

}


new_tree<-tree
droping<-numeric()
if(length(same)>0){
 for(k in 1:length(same)){
   test<-tree$edge[,1]==same[k]
   a<-max(test*tree$edge.length)
   tree$edge.length[test]<-a
   droping[k]<-tree$edge[test,2][1]
 }
   new_tree<-drop.tip(tree, droping)
}


if(length(new_tree$tip.label) == length(tree$tip.label)){break} #ending the do while loop
else{tree <-new_tree}
}   # end of do while


theta<-sum(tree$edge.length) / sum(original_tree$edge.length)
zeta<-(dim(tree$edge)[1]+1)/(dim(original_tree$edge)[1]+1)

output <- data.frame(theta, zeta)
return(output)
}


#######################################################################
### main randomization function for calculating ses and p values
#######################################################################

#'Standardized effect sizes for theta and zeta indices for phylogenetic signal
#'
#'@param cat A dataframe containing data for the categorical trait. The data preferably should be in different columns for each class of the trait and 0 indicating the species does not have this class, while 1 that it does. A species could have more than one categories.
#'@param tree The phylogenetic tree in phylo format. The tree must include all species in the trait dataset
#'@param matching The similarity criterion to be used "complete" for assuming two species as identical if they match in all classes, and "partial" if we assume that species are considered similar if they share at least one class of the trait
#'@param nrep Number of randomizations for estimating the distribution of theta and zeta under random expectations
#'@return The observed values of theta and zeta indices. The values are between 0 and 1, and the smaller the value the stronger the signal. Also by comparison to null model it estimates the Standardized Effect Sizes and the significance level for how likely are the observed values compared to random expectations.
#'
#'@examples
#'tree <- rtree(50)
#'cat <- matrix(sample(c(0,1), 200, replace=TRUE), ncol=4, nrow=50)
#'cat <-data.frame(cat)
#'rownames(cat)<-tree$tip.label
#'ses_theta(cat, tree, "complete", 500)
#'ses_theta(cat, tree, "partial", 500)

ses_theta<-function(cat, tree, matching="complete", nrep){

if(matching!="partial" & matching!="complete"){
stop("Matching criterion for trait comparison among species could be only partial or complete")
}

if(is.data.frame(cat)==FALSE){
stop("Trait data need to be in data frame format, with row names being the species names")
}

#remove lines from trait dataset that have NA
cat<-na.omit(cat)

if(sum(row.names(cat) %in% tree$tip.label)/dim(cat)[1] < 1){
stop("The list of species in trait dataset is not present in tree")
}

if(sum(row.names(cat) %in% tree$tip.label)/length(tree$tip.label) < 1){
tree<-drop.tip(tree,tree$tip.label[-match(row.names(cat), tree$tip.label)])
}

Specs<-row.names(cat)
Traits<-cat

obs<-theta2(Specs,Traits,tree, matching)

dnames<-do.call("cbind", replicate(nrep, Specs, simplify = FALSE))
dnames<-apply(dnames,2,sample)

res<-apply(dnames,2,theta2,Traits,tree,matching)

res<-matrix(unlist(res),nrow=nrep, ncol=2,byrow=TRUE)
means<-apply(res,2,mean)
stds<-apply(res,2,sd)
ses<-(obs-means)/stds
f<-ecdf(res[,1])
p.theta<-f(obs$theta)
f2<-ecdf(res[,2])
p.zeta<-f2(obs$zeta)

outcome<-data.frame(obs$theta, ses$theta, p.theta, obs$zeta, ses$zeta, p.zeta)
return(outcome)
}



