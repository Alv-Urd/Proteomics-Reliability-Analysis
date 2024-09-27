assessment_of_NA_by_day <- function(Dataset_in){ # GitHub
  # NA = Not Assigned value (missing value) in LC-MS/MS data
  # Dataset_in = matrix with proteomics data (including values and NAs)
  # This function returns a matrix: 
    # for each protein, calculates the number of (replicates that are) NAs on any given day
  
  # Assesment of columns in Dataset_in that correspond to biological replicates for each day
  Day_0 = c(1,2,3) # D0R1 was discarded in the original data
  Day_1 = c(4:7)
  Day_2 = c(8:11)
  Day_3 = c(12:15)
  Day_4 = c(16:19)
  Day_5 = c(20:23)
  
  Days = list(Day_0 = Day_0, Day_1 = Day_1, Day_2 = Day_2, Day_3 = Day_3, Day_4 = Day_4, Day_5 = Day_5)
  num_days=6 # Number of time-points (number of columns in the final matrix)
  
  # Matrix: 
    ## as many rows as proteins (rows) in Dataset_in
    ## as many columns as time-points (days, num_days)
  num_NA_by_day = matrix(nrow=dim(Dataset_in)[1], ncol = num_days) 
  rownames(num_NA_by_day) = rownames(Dataset_in)
  All_days = c("Day_0", "Day_1", "Day_2", "Day_3", "Day_4", "Day_5")
  colnames(num_NA_by_day) = All_days
  
  for(prot_index in 1:dim(Dataset_in)[1]){ 
    for(day_index in 1:num_days){
      num_NA_by_day[prot_index,day_index] = sum(is.na(Dataset_in[ prot_index, Days[[day_index]] ])) 
    }
  }
  
  return(num_NA_by_day)
  
}