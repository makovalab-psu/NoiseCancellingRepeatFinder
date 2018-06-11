#===
# see ncrf_positional_filter.py to see how this is used
#===

require(pwr)

do_mx_significance_tests <- function(n,mxFlat,testErrorCounts,effectSize=0.3,power=0.8)
	{ 
	rowLen = length(mxFlat)/n
	for (row in 1:n)
		{
		if (testErrorCounts)
			{
			# this setting of s,e extracts the second half of each row, the error counts
			e = row*rowLen
			s = e-(rowLen-1) + (rowLen/2)
			}
		else
			{
			# this setting of s,e extracts the first half of each row, the match counts
			s = (row-1)*rowLen + 1
			e = s + (rowLen/2) - 1
			}
		cat(assess_significance(mxFlat[s:e],effectSize=effectSize,power=power))
		cat('\n')
		}
	}


assess_significance <- function(case,effectSize=0.3,power=0.8,verbose=F)
	{
	# returns:
	#   TRUE  ==> don't reject null hypothesis; the counts are uniform
	#   FALSE ==> reject null hypothesis;       the counts are biased
	#   NA    ==> unable to run the test

	warning(case)
	if (any(case<5)) return (NA)
	stopifnot(case>=5)
	df <- length(case)-1
	sampleSize <- sum(case)

	alpha <- 0.0000000001 # default alpha value
	tryCatch(
		# Calculating alpha value for fixed power
		alpha <- pwr.chisq.test(w=effectSize,df=df,N=sampleSize,power=power,sig.level=NULL)$sig.level
		,
		# But if an error occurs, we will stick to the default alpha value. 
		# This happens when the alpha value becomes extremely small.
		error=function(error_message)
			{
			message("Uniroot problem: alpha will remain at the minimum value of 0.0000000001.")
			message(error_message)
		  	}
		)

	pvalue <- chisq.test(case)$p.value
	if (verbose)
		write(paste("alpha:",alpha,"pvalue:",pvalue,"df:",df,"sampleSize:",sampleSize), stderr())
	# (pvalue <= alpha) ==> reject null hypothesis;       the counts are biased
	# (pvalue >  alpha) ==> don't reject null hypothesis; the counts are uniform
	return (pvalue>alpha)
	}

