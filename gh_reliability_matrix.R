num_prots = dim(total.log.norm.23.no.NA.R1D0)[1]
num_days=6
##print(num_prots)
##print(num_days)
warning("The total prot num is now redefined after removing those only present in R1D0")
Reliability_matrix = matrix(nrow = num_prots, ncol = num_days)

#****************************************************
# Identify reliable Days 0 - 5
#****************************************************

#Day is reliable if NA <=1


# ALL days
for(prot_index in 1:num_prots){
  #print(paste("protein",prot_index,"---------------------"))
  
  # vector of logical values, saying whether or not the NA_by_day low enough to be  considered a good neighbour:  NA_by_day>=2
  good_neighbour = c(NA_by_day[prot_index,1]<=2, NA_by_day[prot_index,2]<=2, NA_by_day[prot_index,3]<=2, NA_by_day[prot_index,4]<=2, NA_by_day[prot_index,5]<=2, NA_by_day[prot_index,6]<=2)
  # it's a function of day_index! so, it will not be used if all days are reliable on their own, i.e. RD
  
  # loop over all 6 days
  for(day_index in 1:num_days){ # (num_days-1) MUST be in brackets! otherwise R also subtracts 1 from the initial value of i!!!!
    
    # broadly divide into 3 groups, the if statement, according to how many NAs there are in a day:
    
    # - Reliably Detected ----------------------------------------   
    # if up to 1 replicate is an NA  
    if (NA_by_day[prot_index, day_index] <= 1){Reliability_matrix[prot_index,day_index] = 3} # {RD}
    #........................................................
    # if >=2 (2-4) replicates are NAs
    else {
      # if 4 replicates are NAs 
      if (NA_by_day[prot_index, day_index] == 4) {
        # if days are D1-D4, not D5   (D0 cannot be here by definition, as it only has 3 replicates!!!!)
        if (day_index!=6) {  # 
          # if both neighbours are reliable
          if (good_neighbour[day_index-1] & good_neighbour[day_index+1]) {Reliability_matrix[prot_index,day_index] = 1}  # {UU} - NA
          # if it has at least 1 UNRELIABLE neighbour (one or both)! 
          else {Reliability_matrix[prot_index,day_index] = 0}  # {RU} - min_detection  
        }
        # can only be D5!!!! if all 4 replicates are NAs in D5, because it only has one neighbour anyway????
        else {Reliability_matrix[prot_index,day_index] = 0}  # {RU} - min_detection 
      }
      #........................................................
      # if 2-3 replicates are NAs (in any day? i.e. 1-2 replicates are present, except that in the case of D0 this means 0-1 replicates are present!) : RU possible for D0 if 3 NAs; for all days possible UU or UD: 
      else {   
        # D0 when it has all 3 NAs : RU case
        #print(day_index)
        #print(NA_by_day[prot_index,day_index-1])
        if (day_index==1 & NA_by_day[prot_index,day_index]==3) {Reliability_matrix[prot_index,day_index] = 0}  # {RU} - min_detection      
        # for all days, when they have 2-3 NAs: options now are UD or UU
        else {
          # day has at least one good neighbour?
          if (day_index==1) {crit_good_neighbours=good_neighbour[day_index+1]}  # D0
          if (day_index==6) {crit_good_neighbours=good_neighbour[day_index-1]}  # D5
          if (day_index>=2 & day_index<=5) {crit_good_neighbours = (good_neighbour[day_index-1] | good_neighbour[day_index+1])} # D1-4
          
          # Unreliably Detected --------------------------------------------------
          if (crit_good_neighbours) {Reliability_matrix[prot_index,day_index] = 2}    # {UD}    }
          # Unreliably Undetected --------------------------------------------------
          # Day not reliable and neighbours are not reliable
          else {Reliability_matrix[prot_index,day_index] = 1}       #  {UU} - NA   }
        }
      }
    }
    
    
    # Reliably Undetected -------------------------------------------------------
    # CHANGE introduced! to avoid dips:
    
    # when a day has all replicates as NAs (3 in D1 or 4 in all the rest) it is:
    # - Reliably Undetected (min_detection): ONLY if ONE of the 2 neighbours as Unreliable; 
    # - Unreliably Undetected (NA): if BOTH neighbours are Reliable
    # old:
    # Day has all as NAs
    # define the criterion, depending on the day:
    
    # ?????  D1 and D5 have only 1 neighbour - leave them out? apply this only to D1-4
  }
}



#print (Reliability_matrix)
# save("Reliability_matrix", file="Reliability_matrix.RData")
##View(Reliability_matrix)
##summary(Reliability_matrix)
##head(Reliability_matrix)
##colSums(which(Reliability_matrix==0))   - gives Error in colSums(which(Reliability_matrix == 0)) :   'x' must be an array of at least two dimensions
