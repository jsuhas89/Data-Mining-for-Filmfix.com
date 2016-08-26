data<-read.csv("data.csv")
data_cut<-data[2:9]
result<-agnes(t(data_cut),method = "average")
plot(result)