output_header <- function(N, ord, niknots, Se, Sp, beta, delta){
	nbeta <- length(beta)
	sep_lines <- paste0(paste0(rep('-',20+nbeta*7+(nbeta-1)),collapse =''))
	str1 <- paste0(paste0('N                 : '), paste(round(N,4),collapse=" "))
	str2 <- paste0(paste0('ord               : '), paste(round(ord,4),collapse=" "))
	str3 <- paste0(paste0('niknots           : '), paste(round(niknots,4),collapse=" "))
	str4 <- paste0(paste0('Se                : '), paste(round(Se,4),collapse=" "))
	str5 <- paste0(paste0('Sp                : '), paste(round(Sp,4),collapse=" "))
	str6 <- paste0(paste0('true beta         : '), paste(round(beta,4),collapse=" "))
	str7 <- paste0(paste0('true delta        : '), paste(round(delta,4),collapse=" "))
	return(c(sep_lines,str1, str2, str3, str4, str5, str6, str7))
}