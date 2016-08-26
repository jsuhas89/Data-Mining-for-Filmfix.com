#import actual data
data1<-read.csv("data.csv")
data<-data1[2:9]

#initialization
sum<-0
no_of_attributes<-length(colnames(data))
no_of_centroids<-2
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
min<-list()
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
random<-sample(1:length(rownames(data)), no_of_centroids)
for(i in 1:no_of_centroids){
  centroids_initial[[i]]<-data[random[i],]
}


#pick k random centroids from our data set again(these are our actual initializations)
random_again<-sample(1:length(rownames(data)),no_of_centroids)
for(i in 1:no_of_centroids){
  centroids[[i]]<-data[random_again[i],]
}


for(i in 1:no_of_iterations){
  #centroids_list contains k centroids each of which is initialized to 0
  for(m in 1:no_of_centroids){
    centroids_list[[m]]<-matrix(c(0),length(rownames(data)),no_of_attributes)
  }
  
  #initialize the index of the centroids(In R, index starts from 1)
  centroids_index<-rep(c(1),no_of_centroids) 
  
  
  #assign our data points to the relevant centroid
  for(i in 1:length(rownames(data))){
    for(j in 1:no_of_centroids){
      min[j]<-distance(data[i,],centroids[[j]])
    }
    centroids_finder<-which.min(min)
    
    
    #get the ith data into a vector
    for(k in 1:no_of_attributes){    
      idata[k]<-data[i,k] 
    }
    
    centroids_list[[centroids_finder]][centroids_index[centroids_finder],]<-idata
    
    #increment the index to store the next dataset if any    
    centroids_index[centroids_finder]<-centroids_index[centroids_finder] + 1 
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