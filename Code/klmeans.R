#import actual data
bdata_raw<-na.omit(read.csv("meddata_na.csv"))
bdata<-bdata_raw[!duplicated(bdata_raw),][2:10]
adata<-bdata_raw[!duplicated(bdata_raw),]

#initialization
sum<-0
no_of_attributes<-length(colnames(bdata))
no_of_centroids<-4
no_of_iterations<-10
threshold<-0.1
random<-NULL
centroids_initial<-list()
random_again<-NULL
centroids<-list()
centroids_list<-list()
centroids_actual<-list()
centroids_finder<-0
idata<-NULL
jdata<-NULL
centroids_before<-list()
min<-0
mal<-0
ben<-0
TP<-0
FP<-0
PPV<-0

#Euclidean distance function
distance = function(item1,item2){
  return(sqrt(sum((item1 - item2) ^ 2)))
}


#pick initial k random centroids
random<-sample(1:length(rownames(bdata)), no_of_centroids)
for(i in 1:no_of_centroids){
  centroids_initial[[i]]<-bdata[random[i],]
}


#pick k random centroids from our data set again(these are our actual initializations)
random_again<-sample(1:length(rownames(bdata)),no_of_centroids)
for(i in 1:no_of_centroids){
  centroids[[i]]<-bdata[random_again[i],]
}


for(i in 1:no_of_iterations){
  #centroids_list contains k centroids each of which is initialized to 0
  for(m in 1:no_of_centroids){
    centroids_list[[m]]<-matrix(c(0),length(rownames(bdata)),no_of_attributes)
  }
  for(n in 1:no_of_centroids){
    centroids_actual[[n]]<-matrix(c(0),length(rownames(adata)),11)
  }
  
  #initialize the index of the centroids(In R, index starts from 1)
  centroids_index<-rep(c(1),no_of_centroids) 
  
  
  #assign our data points to the relevant centroid
  for(i in 1:length(rownames(bdata))){
    for(j in 1:no_of_centroids){
      min[j]<-distance(bdata[i,],centroids[[j]])
    }
    
    centroids_finder<-which.min(min)
    centroids_finder1<-sort.list(min)
    centroids_finder2<-centroids_finder1[2]
    
    #get the ith data into a vector
    for(k in 1:no_of_attributes){    
      idata[k]<-bdata[i,k] 
    }
    for(l in 1:11){    
      jdata[l]<-adata[i,l] 
    }
    
    centroids_final<-sample(c(centroids_finder,centroids_finder2),1)
    
    
    centroids_list[[centroids_final]][centroids_index[centroids_final],]<-idata
    centroids_actual[[centroids_final]][centroids_index[centroids_final],]<-jdata
    
    #same code for second closest centroid
    #centroids_list[[centroids_finder2]][centroids_index[centroids_finder2],]<-idata
    #centroids_actual[[centroids_finder2]][centroids_index[centroids_finder2],]<-jdata
    
    #increment the index to store the next dataset if any    
    centroids_index[centroids_final]<-centroids_index[centroids_final] + 1 
    
    #centroids_index[centroids_finder2]<-centroids_index[centroids_finder2] + 1
  }
  
  #decrement the index by 1 since the previous step would have incremented by 1 during last iteration
  for(i in 1:length(centroids_index)){
    centroids_index[i]<-centroids_index[i]-1
  }
  
  
  #calculate the new set of centroids
  centroids_before<-centroids
  for(i in 1:no_of_centroids){ 
    for(j in 1:no_of_attributes){
      sum<-0
      for(k in 1:(centroids_index[i])){
        sum<-sum + centroids_list[[i]][k,j]
      }
      sum<-sum/(centroids_index[i])
      centroids[[i]][j]<-sum
    }
  }
  
  
  #calculate f0 as per the algorithm
  f0<-0
  for(i in 1:length(centroids_before)){ 
    f0_each<-0
    for(j in 1:length(centroids_initial)){
      f0_each<-f0_each + distance(centroids_before[[i]],centroids_initial[[j]])
    }  
    f0<-f0 + f0_each
  }
  
  
  #calculate f1 as per the algorithm  
  f1<-0
  for(i in 1:length(centroids)){
    f1_each<-0
    for(j in 1:length(centroids_before)){
      f1_each<-f1_each + distance(centroids[[i]],centroids_before[[j]])
    }
    f1<-f1 + f1_each
  }
  change_per<-(f1-f0)/f0
  if(abs(change_per)<=threshold){ 
    break
  }
  centroids_initial<-centroids_before 
}

#print centroid indexes
print(centroids_index)


#code to calculate the total PPV of the clusters
for(k in 1:no_of_centroids){
  mal<-length(which(centroids_actual[[k]][,11] == 4))
  print(mal)
  ben<-length(which(centroids_actual[[k]][,11] == 2))
  print(ben)
  if(ben>=mal)
  { 
    TP = TP + ben
    FP = FP + mal
    #PPV = TP/(TP+FP)
    #cat(sprintf("PPV of centroid %f is: %f\n",k,PPV))
  }
  else
  {
    TP = TP + mal
    FP = FP + ben
    #PPV = TP/(TP+FP)
    #cat(sprintf("PPV of centroid %f is: %f\n",k,PPV))
  }
}
PPV = TP/(TP+FP)
cat(sprintf("PPV is: %f\n",PPV))