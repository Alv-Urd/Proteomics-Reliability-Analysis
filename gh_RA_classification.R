RA_classification = function(Dataset_in){ # GitHub
  # Dataset_in = NA matrix 
  
 CLASSIFICATION <-  Dataset_in %>%
    mutate(D0 = ifelse(Day_0 <= 1, "RD",
                       ifelse(Day_0 == 3, "RU",
                              ifelse(Day_0 >= 2 & Day_1 <= 2, "UD",
                                     ifelse(Day_0 >= 2 & Day_1 > 2, "UU",
                                            "Error"
                                     )))),
           D1 = ifelse(Day_1 <= 1, "RD",
                       ifelse(Day_1 == 4 & (Day_0 >= 2 | Day_2 >= 2), "RU",
                              ifelse(Day_1 >= 2 & (Day_0 <= 2 | Day_2 <= 2), "UD",
                                     ifelse(Day_1 >= 2 & (Day_0 > 2 | Day_2 > 2), "UU",
                                            "Error"
                                     )))),
           D2 = ifelse(Day_2 <= 1, "RD",
                       ifelse(Day_2 == 4 & (Day_1 >= 2 | Day_3 >= 2), "RU",
                              ifelse(Day_2 >= 2 & (Day_1 <= 2 | Day_3 <= 2), "UD",
                                     ifelse(Day_2 >= 2 & (Day_1 > 2 | Day_3 > 2), "UU",
                                            "Error"
                                     )))),
           D3 = ifelse(Day_3 <= 1, "RD",
                       ifelse(Day_3 == 4 & (Day_2 >= 2 | Day_4 >= 2), "RU",
                              ifelse(Day_3 >= 2 & (Day_2 <= 2 | Day_4 <= 2), "UD",
                                     ifelse(Day_3 >= 2 & (Day_2 > 2 | Day_4 > 2), "UU",
                                            "Error"
                                     )))),
           D4 = ifelse(Day_4 <= 1, "RD",
                       ifelse(Day_4 == 4 & (Day_3 >= 2 | Day_5 >= 2), "RU",
                              ifelse(Day_4 >= 2 & (Day_3 <= 2 | Day_5 <= 2), "UD",
                                     ifelse(Day_4 >= 2 & (Day_3 > 2 | Day_5 > 2), "UU",
                                            "Error"
                                     )))), 
           D5 = ifelse(Day_5 <= 1, "RD",
                       ifelse(Day_5 == 4, "RU",
                              ifelse(Day_5 >= 2 & Day_4 <= 2, "UD",
                                     ifelse(Day_5 >= 2 & Day_4 > 2, "UU",
                                            "Error"
                                     ))))
    )
 Dataset_out <- CLASSIFICATION[,7:13]
 return(Dataset_out)
}