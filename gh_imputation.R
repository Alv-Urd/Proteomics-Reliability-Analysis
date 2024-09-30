# Original proteomics data comes from the following file: total.log.norm.23.no.NA.R1D0.Reliability.no.undetected.csv

prot <- read.table(paste0(input_data,"total.log.norm.23.no.NA.R1D0.Reliability.no.undetected.csv"), sep=",", header=TRUE)


## Lists of proteins with leading or trailing NAs (i.e.: all replicates NA at time 0, or times 0 & 1, or times 0 & 1 & 2,..., or times 4 & 5, or time 5)

#a) Counting leading NAs:

t0na<-which(is.na(prot$R2D0) & is.na(prot$R3D0) & is.na(prot$R4D0))
#Proteins with NA in all replicates at time 0: 498

t1na<-which(is.na(prot$R2D0) & is.na(prot$R3D0) & is.na(prot$R4D0) & is.na(prot$R1D1) & is.na(prot$R2D1) & is.na(prot$R3D1) & is.na(prot$R4D1))
# , of which 47 have also NAs in all replicates at time 1

t2na<-which(is.na(prot$R2D0) & is.na(prot$R3D0) & is.na(prot$R4D0) & is.na(prot$R1D1) & is.na(prot$R2D1) & is.na(prot$R3D1) & is.na(prot$R4D1) & is.na(prot$R1D2) & is.na(prot$R2D2) & is.na(prot$R3D2) & is.na(prot$R4D2))
# , of which 13 have NAs in all replicates at time 2. There are no cases with all replicates at times 0, 1, 2, and 3

#b) Counting trailing NAs:

t5na<-which(is.na(prot$R1D5) & is.na(prot$R2D5) & is.na(prot$R3D5) & is.na(prot$R4D5))
# Proteins with NA in all replicates at time 0: 678

t4na<-which(is.na(prot$R1D5) & is.na(prot$R2D5) & is.na(prot$R3D5) & is.na(prot$R4D5) & is.na(prot$R1D4) & is.na(prot$R2D4) & is.na(prot$R3D4) & is.na(prot$R4D4))
# , of which 156 have also NAs in all replicates at time 4

t543na<-which(is.na(prot$R1D5) & is.na(prot$R2D5) & is.na(prot$R3D5) & is.na(prot$R4D5) & is.na(prot$R1D4) & is.na(prot$R2D4) & is.na(prot$R3D4) & is.na(prot$R4D4) & is.na(prot$R1D3) & is.na(prot$R2D3) & is.na(prot$R3D3) & is.na(prot$R4D3))
# , of which 47 have also NAs in all replicates at time 3

t5432na<-which(is.na(prot$R1D5) & is.na(prot$R2D5) & is.na(prot$R3D5) & is.na(prot$R4D5) & is.na(prot$R1D4) & is.na(prot$R2D4) & is.na(prot$R3D4) & is.na(prot$R4D4) & is.na(prot$R1D3) & is.na(prot$R2D3) & is.na(prot$R3D3) & is.na(prot$R4D3)& is.na(prot$R1D2) & is.na(prot$R2D2) & is.na(prot$R3D2) & is.na(prot$R4D2))
# , of which 15 have also NAs in all replicates at time 2. There are no cases with NAs in all replicates at times 1, 2, 3, 4, 5.


###NAs Imputation

t05na<-sort(union(t0na,t5na))
# 'union' discard repeated values. t05na are those proteins with 1 NA at time 0 or 1 NA at time 5 (or both).

# Computing median values per time point and protein:
t0<-prot[,c(1,2,3,4)]
rownames(t0)<-t0[,1]
t0<-t0[,2:4]

t0[,4] <- rowMedians(as.matrix(t0), na.rm=TRUE)

t1<-prot[,c(1,5,6,7,8)]
rownames(t1)<-t1[,1]
t1<-t1[,2:5]

t1[,5] <- rowMedians(as.matrix(t1), na.rm=TRUE)


t2<-prot[,c(1,9,10,11,12)]
rownames(t2)<-t2[,1]
t2<-t2[,2:5]

t2[,5] <- rowMedians(as.matrix(t2), na.rm=TRUE)


t3<-prot[,c(1,13,14,15,16)]
rownames(t3)<-t3[,1]
t3<-t3[,2:5]

t3[,5] <- rowMedians(as.matrix(t3), na.rm=TRUE)


t4<-prot[,c(1,17,18,19,20)]
rownames(t4)<-t4[,1]
t4<-t4[,2:5]

t4[,5] <- rowMedians(as.matrix(t4), na.rm=TRUE)

t5<-prot[,c(1,21,22,23,24)]
rownames(t5)<-t5[,1]
t5<-t5[,2:5]

t5[,5] <- rowMedians(as.matrix(t5), na.rm=TRUE)


au <- cbind(t0[,"V4"], t1[,"V5"], t2[,"V5"], t3[,"V5"], t4[,"V5"], t5[,"V5"])
rownames(au)<-rownames(t5)
colnames(au)<-c("0","1","2","3","4","5")
# 'au' contains as many columns as data points.
# The data is the medians of the replicates.
# But a day without no data produces NaN. These are the points we need to spline-interpolate. Or else, knn-interpolate.

# Replacing NaN for NA:
au[is.nan(au)] <- NA

sum(!complete.cases(au))
# 1,541 proteins with at least one time point without data.

## Second approach: knn-imputation:

# The following is to obtain the matrix with NA in non-empty days imputed as median for the day, while empty days remains empty:

t0cc<-t0
t1cc<-t1
t2cc<-t2
t3cc<-t3
t4cc<-t4
t5cc<-t5

t0cc$R2D0= coalesce(t0$R2D0, au[,"0"])
t0cc$R3D0= coalesce(t0$R3D0, au[,"0"])
t0cc$R4D0= coalesce(t0$R4D0, au[,"0"])

t1cc$R1D1= coalesce(t1$R1D1, au[,"1"])
t1cc$R2D1= coalesce(t1$R2D1, au[,"1"])
t1cc$R3D1= coalesce(t1$R3D1, au[,"1"])
t1cc$R4D1= coalesce(t1$R4D1, au[,"1"])

t2cc$R1D2= coalesce(t2$R1D2, au[,"2"])
t2cc$R2D2= coalesce(t2$R2D2, au[,"2"])
t2cc$R3D2= coalesce(t2$R3D2, au[,"2"])
t2cc$R4D2= coalesce(t2$R4D2, au[,"2"])

t3cc$R1D3= coalesce(t3$R1D3, au[,"3"])
t3cc$R2D3= coalesce(t3$R2D3, au[,"3"])
t3cc$R3D3= coalesce(t3$R3D3, au[,"3"])
t3cc$R4D3= coalesce(t3$R4D3, au[,"3"])

t4cc$R1D4= coalesce(t4$R1D4, au[,"4"])
t4cc$R2D4= coalesce(t4$R2D4, au[,"4"])
t4cc$R3D4= coalesce(t4$R3D4, au[,"4"])
t4cc$R4D4= coalesce(t4$R4D4, au[,"4"])

t5cc$R1D5= coalesce(t5$R1D5, au[,"5"])
t5cc$R2D5= coalesce(t5$R2D5, au[,"5"])
t5cc$R3D5= coalesce(t5$R3D5, au[,"5"])
t5cc$R4D5= coalesce(t5$R4D5, au[,"5"])

prot_3 <- cbind (t0cc[,-4],t1cc[,-5],t2cc[,-5],t3cc[,-5],t4cc[,-5],t5cc[,-5])
# So, 'prot_3' replace Nas in time points with at least one valid replicate, with the median for that time point.
# But leaves empty time points (i.e. all replicates NA), as NA.


#knn-imputation:

library(DMwR2)


prot_3_knn<-knnImputation(prot_3)
# 'Prot_3_knn' knn-imputes empty days and median-impute NAs in non-empty days.


# Choice of k is very critical – A small value of k means that noise will have a higher influence on the result. 
# A large value make it computationally expensive and kinda defeats the basic philosophy behind KNN
# (that points that are near might have similar densities or classes ).
# A simple approach to select k is set k = n^(1/2). In my case is 7. Many sources say it needs to be an odd number (?).
# I leave K =10, which is the default.

anyNA(prot_3_knn)

write.table(prot_3_knn, "prot_3_knn", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)


#Now, make sure no value is below 15 (minimum value):
prot_3_knn[rowSums(prot_3_knn < 14.99999) != 0,] 

