library(AlphaSimR)

## Initial parameters
nt     <-  5 # number of traits
ngener <- 10 # number of simulation generations/years


varG <- c(27.60, 37.00, 48.16, 9.00, 847.60)
resVar <- matrix(c(32.40,    12.12,     9.46,     3.21,    31.85,
                   12.12,    37.00,    15.35,     2.84,    42.54,
		   9.46,    15.35,    37.84,     3.47,    58.51,
		   3.21,     2.84,     3.47,     6.00,    11.65,
		   31.85,    42.54,    58.51,    11.65,   782.40), ncol = nt)   # residual variance covariance matrix


genCor <- matrix(c(1.0,       0.970077,  0.950124,  0.653524,  0.770839,
                   0.970077,  1.0,       0.980035,  0.701436,  0.84002,
		   0.950124,  0.980035,  1.0,       0.680141,  0.850324,
		   0.653524,  0.701436,  0.680141,  1.0,       0.59079,
		   0.770839,  0.84002,   0.850324,  0.59079,   1.0), ncol = nt) # Genetic correlation matrix
 
gVar <- matrix(c(27.60, 31.00, 34.64, 10.30, 117.90,
                 31.00, 37.00, 41.37, 12.80, 148.76,
                 34.64, 41.37, 48.16, 14.16, 171.80,
		 10.30, 12.80, 14.16, 9.00,  51.60,
		 117.90, 148.76, 171.80, 51.60, 847.60), ncol = nt) # genetic variance covariance matrix



w  = c(1, 1, 1, 1, 1)     # Trait economic weights
cat("Eco weights=",w,"\n")


nsel_female <- 1000       # number of selected females
nsel_male   <- 100        # number of selected males
nbase       <- 2*nsel_female # base population size, has about equal number of both females and males
Ne          <- 200        # effective population size 

cat("\nBase size       :", nbase,"\n")
cat("Selected males  :", nsel_male,"\n")
cat("Selected females:", nsel_female,"\n")

Nchr   <- 30              # number of chromosomes
nsites <- 3600            # number of segregation sites
jmp    <- 2               # SNP jump to segregation sites
Qjmp   <- 120             # QTL jump to segregation sites
Nsnp   <- Nchr*nsites/jmp # number of SNP markers
NsnpPerChr <- Nsnp/Nchr   # number of SNP markers per chromosome
NQTL       <- Nchr*nsites/Qjmp # number of QTLs
NQTLperChr <- max(NQTL/Nchr,1) # number of QTLs per chromosome

cat("\nN chromosomes   :", Nchr,"\n")
cat("N segreg. sites   :", nsites,"\n")
cat("N SNP             :", Nsnp,"\n")
cat("N SNP per chr     :", NsnpPerChr,"\n")
cat("N QTL             :", NQTL,"\n")
cat("N QTL per chr     :", NQTLperChr,"\n")

# Create founder haplotypes
tmp = "founderPopRep.RData"       # haplotypes folder
if (!file.exists(tmp)) {                                    # don't create new haplotypes if they already exist
  cat("\nrunMacs...\n")
  founderPop = runMacs(nInd = Ne,        # Number of individuals
                     nChr = Nchr,        # Number of chromosomes
                     segSites = nsites,  # Number of segregating sites to keep per chromosome
                     species = "CATTLE") # The CATTLE history
  save(founderPop, file=tmp)
} else {
  cat("\nload founder population...\n")
  load(tmp)
  print(tmp)
}

# If simulating more than one trait, all traits will be pleiotrophic with correlated additive effects
cat("SimParam...\n")
# Set simulation parameters
SP = SimParam$
  new(founderPop)$
  addTraitA(NQTLperChr, mean=c(0,0,0,0,0), var=varG, corA=genCor,  # number of QTL per chr; genetic var&correlation
            gamma=TRUE, shape=0.4)$ # gamma density with shape leads to large QTLs
  setVarE(varE=resVar)$    # residual covariance matrix
  setSexes("yes_sys")$    # randomly assign a sex to each individual
  addSnpChip(NsnpPerChr, minSnpFreq=0.05)$  # Randomly assigns eligible SNPs to a SNP chip
  setTrackPed(TRUE)        # pedigree tracking ON
 
cat("newPop...\n")
pop = newPop(founderPop, simParam=SP) # Create new population
 
pop_base = randCross(pop = pop, nCrosses = nbase)  #Generating 6000 individuals

pop_all <- c(pop, pop_base) 
pop     <- pop_base         
cat("Pop so far, N indiv(pop_all)=", nInd(pop_all),"\n")
cat("Breeding pop, N indiv(pop)  =", nInd(pop),"\n")

gmean   <- meanG(pop) # Returns the mean genetic values for all traits
gMean_1 <- gmean[1]
gMean_2 <- gmean[2]
gMean_3 <- gmean[3]
gMean_4 <- gmean[4]
gMean_5 <- gmean[5]
cat("The mean genetic values is ", gmean, "\n")
 
id_now <- as.numeric(slot(pop_all,"id"))
nn     <- length(id_now)
gener  <- matrix(1,nn,1)
id_num <- matrix(id_now, nn,1)


for(year in 2:(ngener+1)){
  cat("Year...", year,"N ind=",nInd(pop_all),"\n")
  pop_cross   <- randCross(pop, nCrosses=nsel_female) # Make random crosses
  id_now      <- as.numeric(slot(pop_cross,"id"))
  nn          <- length(id_now)
  taulu       <- matrix(year,nn,1)
  gener       <- c(gener, taulu)
  id_num      <- c(id_num, id_now)
  gmean       <- meanG(pop_cross)                     # Mean genetics of the crosses
  # sukup       <- slot(pop_cross, "sex")
  gMean_1     = c(gMean_1, gmean[1])                  # Add the mean by trait
  gMean_2     = c(gMean_2, gmean[2])
  gMean_3     = c(gMean_3, gmean[3])  
  gMean_4     = c(gMean_4, gmean[4])
  gMean_5     = c(gMean_5, gmean[5])  
  pop_tot     <- c(pop, pop_cross)                    # current population and new crosses
  pop_all     <- c(pop_all, pop_cross)                # complete population
  pop_males   <- selectInd(pop_tot, nInd=nsel_male,   sex="M", trait=selIndex, b=w) # male select
  pop_females <- selectInd(pop_tot, nInd=nsel_female, sex="F", trait=selIndex, b=w) # female select
  pop         <- c(pop_males, pop_females)            # new current population
}

cat("Final     , N indiv(pop_all)=", nInd(pop_all),"\n")
cat("Writing files...\n")
sp             <- slot(pop_all, "sex") # get sex of individuals in the pedigree
n              <- length(sp)
sukup          <- matrix(1,n,1)
sukup[sp=="F"] <- 2

pedi    <- as.data.frame(cbind(getPed(pop_all), gener, sukup))# pedigree of population individuals
id      <- as.matrix(pedi$id, n,1)
bvs     <- as.data.frame(cbind(id, bv(pop_all)))     # breeding values of individuals
geno    <- as.data.frame(cbind(id,pullSnpGeno(pop_all, simParam = SP))) # SNP genotype data
pheno   <- as.data.frame(cbind(id,gener, sukup, pheno(pop_all)))   # phenotypes
QTLgeno <- as.data.frame(cbind(id,pullQtlGeno(pop_all, simParam = SP))) # QTLs

write.table(pedi, file="sim.ped", quote=FALSE, col.names = FALSE, row.names = FALSE)
write.table(bvs, file="TBV.txt", quote=FALSE, col.names = FALSE, row.names = FALSE)
write.table(pheno, file="pheno.txt", quote=FALSE, col.names = FALSE, row.names = FALSE)
write.table(QTLgeno, file="QTL.txt", quote=FALSE, col.names = FALSE, row.names = FALSE)

# get the QTL values
h      <- SP$traits
ade1   <- slot(h[[1]],"addEff")
ade2   <- slot(h[[2]],"addEff")
ade3   <- slot(h[[3]],"addEff")
ade4   <- slot(h[[4]],"addEff")
ade5   <- slot(h[[5]],"addEff")
QTLeff <- cbind(ade1, ade2, ade3, ade4, ade5)
write.table(QTLeff, file="QTLeff.txt", quote=FALSE, col.names = FALSE, row.names = FALSE)
write.table(geno, file="SNP.txt", quote=FALSE, col.names = FALSE, row.names = FALSE)

cat("Done.\n")





























