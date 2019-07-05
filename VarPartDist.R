Dist_Matrix: A dissimilarity matrix
Env: Environmental table with samples as rows and environmental variables as columns.
Geo: Geographic table with samples as rows and geographic/cartesian coordinates as columns.
#If Geographic table is geographic coordinate, Geo_Co=TRUE.
#If Geographic table is cartesian coordinate, Geo_Co=FALSE.
Number_Permutations: the number of permutations required in Permutation Test

VarPartDist<-function (Dist_Matrix,Env,Geo,Geo_Co=TRUE,Number_Permutations=999) {
  require(SoDA)
  require(geosphere)
  require(vegan)

  Env<-as.data.frame(Env)
  Geo<-as.data.frame(Geo)
  if(Geo_Co==TRUE){
         geo<-geoXY(Geo[,1],Geo[,2])########Translate the Geographic coordinate to Cartesian coordinates
  } else {
         geo<-Geo
  }   
  pc<-pcnm(dist(geo))$vectors
  pc<-as.data.frame(pc)

  ##Test and forward selection of PCNM variables
  mod1<-capscale(Dist_Matrix~.,Env,add=T)
  mod0<-capscale(Dist_Matrix~1,Env,add=T)
  mod<-ordiR2step(mod0,scope=formula(mod1),perm.max=999)
  Env.P<-anova.cca(mod,permutations=Number_Permutations)[[4]][1]
  Env_se<-mod$CCA$biplot
  Env_se<-Env[,rownames(Env_se)]


  ##Test and forward selection of environmental variables
  mod1<-capscale(Dist_Matrix~.,pc,add=T)
  mod0<-capscale(Dist_Matrix~1,pc,add=T)
  mod<-ordiR2step(mod0,scope=formula(mod1),perm.max=999)
  Geo.P<-anova.cca(mod,permutations=Number_Permutations)[[4]][1]
  Geo_se<-mod$CCA$biplot
  Geo_se<-pc[,rownames(Geo_se)]

  ##variation partitioning
  Var_mod<-varpart(Dist_Matrix,Env_se,Geo_se)
  Env.R2<-Var_mod[1]$`part`[[2]][1,3]
  Geo.R2<-Var_mod[1]$`part`[[2]][2,3]
  Spe.R2<-Var_mod[1]$`part`[[3]][1,3]
  Dis.R2<-Var_mod[1]$`part`[[3]][3,3]
  Com.R2<-Var_mod[1]$`part`[[3]][2,3]
  Spe.P<-anova.cca(capscale(Dist_Matrix~.+Condition(as.matrix(Geo_se)),data=Env_se),permutations=Number_Permutations)[[4]][1]
  Dis.P<-anova.cca(capscale(Dist_Matrix~.+Condition(as.matrix(Env_se)),data=Geo_se),permutations=Number_Permutations)[[4]][1]

  result <- mat.or.vec(5,2)
  result[1,1] <- Env.R2
  result[2,1] <- Geo.R2
  result[3,1] <- Spe.R2
  result[4,1] <- Dis.R2
  result[5,1] <- Com.R2
  
  result[1,2] <- Env.P
  result[2,2] <- Geo.P
  result[3,2] <- Spe.P
  result[4,2] <- Dis.P
  result[5,2] <- NA
  colnames(result) <- c("Adj.R.squared","p-value")
  rownames(result) <- c("Enviromental","Geographic","Species sorting","Dispersal limitation","Combinated fraction")
  return(result)
}



