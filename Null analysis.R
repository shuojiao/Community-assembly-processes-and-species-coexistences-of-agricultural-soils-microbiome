#### R codes for null model analysis modified according to Stegen et al. (2013) ####

phylo: Phylogenetic tree of each OTU
comun: A community table with samples as rows and OTUs as columns. 

#Beta_NTI
Beta_NTI<-function(phylo,comun,beta.reps=999){
   require(picante)

   comun=t(comun)
   match.phylo.comun = match.phylo.data(phylo, t(comun))
   beta.mntd.weighted = as.matrix(comdistnt(t(match.phylo.comun$data),cophenetic(match.phylo.comun$phy),abundance.weighted=T))

   rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(match.phylo.comun$data),ncol(match.phylo.comun$data),beta.reps))
   for (rep in 1:beta.reps) {
       rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(match.phylo.comun$data),taxaShuffle(cophenetic(match.phylo.comun$phy)),abundance.weighted=T,exclude.conspecifics = F))
       print(c(date(),rep))
   }

   weighted.bNTI = matrix(c(NA),nrow=ncol(match.phylo.comun$data),ncol=ncol(match.phylo.comun$data))
   for(columns in 1:(ncol(match.phylo.comun$data)-1)) {
          for(rows in (columns+1):ncol(match.phylo.comun$data)) {
                  rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
                  weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals)
                  rm("rand.vals")
          }
   }
   
  rownames(weighted.bNTI) = colnames(match.phylo.comun$data);
  colnames(weighted.bNTI) = colnames(match.phylo.comun$data);
  return(as.dist(weighted.bNTI))
}

#RC_bray
raup_crick= function(comun, reps=999){
  ## count number of sites and total species richness across all plots (gamma)
  n_sites<-nrow(comun)
  gamma<-ncol(comun)
  ##build a site by site matrix for the results, with the names of the sites in the row and col names:
  results<-matrix(data=NA, nrow=n_sites, ncol=n_sites, dimnames=list(row.names(comun), row.names(comun)))
  ##make the comun matrix into a new, pres/abs. matrix:
  ceiling(comun/max(comun))->comun.inc
  ##create an occurrence vector- used to give more weight to widely distributed species in the null model:
  occur<-apply(comun.inc, MARGIN=2, FUN=sum)
  ##create an abundance vector- used to give more weight to abundant species in the second step of the null model:
  abundance<-apply(comun, MARGIN=2, FUN=sum)
  ##make_null:
  ##looping over each pairwise community combination:
  for(null.one in 1:(nrow(comun)-1)){
    for(null.two in (null.one+1):nrow(comun)){
      null_bray_curtis<-NULL
      for(i in 1:reps){
        ##two empty null communities of size gamma:
        com1<-rep(0,gamma)
        com2<-rep(0,gamma)
        ##add observed number of species to com1, weighting by species occurrence frequencies:
        com1[sample(1:gamma, sum(comun.inc[null.one,]), replace=FALSE, prob=occur)]<-1
        com1.samp.sp = sample(which(com1>0),(sum(comun[null.one,])-sum(com1)),replace=TRUE,prob=abundance[which(com1>0)]);
        com1.samp.sp = cbind(com1.samp.sp,1); # head(com1.samp.sp);
        com1.sp.counts = as.data.frame(tapply(com1.samp.sp[,2],com1.samp.sp[,1],FUN=sum)); colnames(com1.sp.counts) = 'counts'; # head(com1.sp.counts);
        com1.sp.counts$sp = as.numeric(rownames(com1.sp.counts)); # head(com1.sp.counts);
        com1[com1.sp.counts$sp] = com1[com1.sp.counts$sp] + com1.sp.counts$counts; # com1;
        #sum(com1) - sum(spXsite[null.one,]); ## this should be zero if everything work properly
        rm('com1.samp.sp','com1.sp.counts');			
        ##same for com2:
        com2[sample(1:gamma, sum(comun.inc[null.two,]), replace=FALSE, prob=occur)]<-1
        com2.samp.sp = sample(which(com2>0),(sum(comun[null.two,])-sum(com2)),replace=TRUE,prob=abundance[which(com2>0)]);
        com2.samp.sp = cbind(com2.samp.sp,1); # head(com2.samp.sp);
        com2.sp.counts = as.data.frame(tapply(com2.samp.sp[,2],com2.samp.sp[,1],FUN=sum)); colnames(com2.sp.counts) = 'counts'; # head(com2.sp.counts);
        com2.sp.counts$sp = as.numeric(rownames(com2.sp.counts)); # head(com2.sp.counts);
        com2[com2.sp.counts$sp] = com2[com2.sp.counts$sp] + com2.sp.counts$counts; # com2;
        # sum(com2) - sum(spXsite[null.two,]); ## this should be zero if everything work properly
        rm('com2.samp.sp','com2.sp.counts');
        null.comun = rbind(com1,com2); # null.comun;
        ##calculate null bray curtis
        null_bray_curtis[i] = distance(null.comun,method='bray-curtis');
      }; # end reps loop
      ## empirically observed bray curtis
      obs.bray = distance(comun[c(null.one,null.two),],method='bray-curtis');
      ##how many null observations is the observed value tied with?
      num_exact_matching_in_null = sum(null_bray_curtis==obs.bray);
      ##how many null values are smaller than the observed *dissimilarity*?
      num_less_than_in_null = sum(null_bray_curtis<obs.bray);
      rc = ((num_less_than_in_null +(num_exact_matching_in_null)/2)/reps)
      ##modification of raup crick standardizes the metric to range from -1 to 1 instead of 0 to 1
      rc = (rc-.5)*2
      results[null.two,null.one] = round(rc,digits=2); ##store the metric in the results matrix
      print(c(null.one,null.two,date()));
    }; ## end null.two loop
  }; ## end null.one loop

  results<-as.dist(results)
  return(results)
}

