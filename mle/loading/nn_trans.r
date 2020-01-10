nn_trans <- function(vec,log_trans=FALSE){
	nvec = length(vec)
	if(log_trans){
		vec[2:(nvec/2)]<-log(vec[2:(nvec/2)])
		vec[(nvec/2+2):nvec]<-log(vec[(nvec/2+2):nvec])
	}else{
		vec[2:(nvec/2)]<-exp(vec[2:(nvec/2)])
		vec[(nvec/2+2):nvec]<-exp(vec[(nvec/2+2):nvec])
	}
	return(vec)
}