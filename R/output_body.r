output_body <- function(costs, err, nest_err, next_beta, next_delta){
	nbeta <- length(next_beta)
	sep_lines <- paste0(paste0(rep('-',20+nbeta*7+(nbeta-1)),collapse =''))
	str1 <- paste0(paste0('last three costs  : '), paste(round(tail(costs,3),3),collapse=" "))
	str2 <- paste0(paste0('costs change rate : '), paste0(round(err,4)*100, '%'))
	str3 <- paste0(paste0('nest cost diff    : '), paste(round(nest_err,4),collapse=" "))
	str4 <- paste0(paste0('est beta          : '), paste(round(next_beta,4),collapse=" "))
	str5 <- paste0(paste0('est delta         : '), paste(round(next_delta,4),collapse=" "))
	return(c(sep_lines, str1, str2, str3, str4, str5))
}