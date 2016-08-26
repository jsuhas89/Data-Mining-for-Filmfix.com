data_knn<-read.csv("data_knn.csv")
data_knn1<-data_knn[2:11]
cl<-matrix(1:14,14,1)
data_test<-read.csv("data_test.csv")
data_test1<-data_test[2:11]
knn(data_knn1,data_test1,cl,k=1)