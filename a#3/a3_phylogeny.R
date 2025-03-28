# BIOL 365 Assign3 - Phylogenetic analysis in R
# < Hafsa Mueen >, 20947674 
R.version.string
date()

# loading packages + setting up work directory 
library(msa)
library(phangorn)
setwd("~/Documents/BIOL365/a#3")

# loading in FASTA sequences for COX1 gene for the given species:
# seals, walrus, sea lion plus outgroups cat and dog 
pin1_data = readDNAStringSet("~/Documents/BIOL365/a#3/Pinnipeds_initdata.txt")

# creating the multiple sequence alignment of the given species
pin1_msa = msa(pin1_data)
print(pin1_msa, show = "complete") # confirming alignment 
pin1_phydat = as.phyDat(pin1_msa)  # correctly formatting for phangorn 

# analyzing how initial species are related to each other; 
# no optimized evolutionary model used 
pin1_dm = dist.ml(pin1_phydat)                  # initial distance matrix 
pin1_nj_tree = NJ(pin1_dm)                      # creating NJ tree
plot.phylo(pin1_nj_tree, cex = 0.4)             # plotting rooted NJ tree
plot.phylo(pin1_nj_tree, "unrooted", cex = 0.4) # plotting unrooted NJ tree 

# creation of the new pinniped MSA with the addition of marine mammals,
# manatees, seacows, two closely related species and outgroup non-mammalians
pin_new_data = readDNAStringSet("~/Documents/BIOL365/a#3/Pinnipeds_initdata_new1.txt") 
pin_new_msa = msa(pin_new_data) 
sink(file = "new_msa_output.txt")       # saving output to text file
print(pin_new_msa, show = "complete")   # confirming alignment 
sink() 
pin_new_phydat = as.phyDat(pin_new_msa) # correctly formatting for phangorn

# analyzing how initial and newly added species are related to each other
pin_new_dm = dist.ml(pin_new_phydat) 
pin_new_nj_tree = NJ(pin_new_dm)
plot.phylo(pin_new_nj_tree, cex = 0.4) # plotting rooted NJ tree

# refining tree by choosing best evolutionary model to best fit data 
mt = modelTest(pin_new_phydat, model = c("JC", "F81", "K80", "HKY", "SYM", "GTR")) # evolutionary models used
mt                # printing all results of the model test 
which.min(mt$BIC) # printed out 24 

# choosing the model with the lowest BIC 
env = attr(mt, "env")
bestmodel = mt$Model[which.min(mt$BIC)]
print(bestmodel)                          # checking which model is the best
fitStart = eval(get(bestmodel, env), env)
fit = optim.pml(fitStart, rearrangement = "stochastic", optGamma=TRUE, optInv=TRUE, bestmodel)
fit

# adding bootstrap values to the now fitted tree and plotting the tree
bs = bootstrap.pml(fit, bs = 100, optNni = TRUE)
print(bs)
plotBS(midpoint(fit$tree), bs, p = 0, type = "p", cex = 0.6)
            

