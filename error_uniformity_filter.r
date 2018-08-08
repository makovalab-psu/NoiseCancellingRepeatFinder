#===
# see error_uniformity_filter.py to see how this is used
#
# Note that this is only used ito support the error_uniformity_filter.py option
# --method=chi-squared. This method has been found to be flawed, and is not
# recommended.
#===

require(pwr)

do_mx_significance_tests <- function(n,mxFlat,testErrorCounts,effectSize=0.3,power=0.8,verbose=F)
	{
	# n
	#	is the number of flattened rows (the number of tests to perform)
	# mxFlat
	#	is a flattened matrix of the test values; each flattened row has 2*m
	#	columns, where m is the vector size of a test; thus each flattened row
	#   has two test vectors, which for our original purposes were "match
	#	counts" and "error counts"
	# testErrorCounts
	#	is TRUE  if we're to test the second test vector in each flattened row
	#	is FALSE if we're to test the first  test vector in each flattened row
	#   in either case, only one of the two vectors is tested

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
		if (verbose)
			write(paste("row",row,"events:",paste(mxFlat[s:e],collapse=",")), stderr())
		result = assess_significance(mxFlat[s:e],effectSize=effectSize,power=power)
		if (verbose) # nb: assess_significance's verbose output fails to get out
			write(paste("alpha:",result$alpha,"pvalue:",result$pvalue,"df:",result$df,"sampleSize:",result$sampleSize), stderr())
		cat(result$outcome)
		cat('\n')
		}
	}


assess_significance <- function(case,effectSize=0.3,power=0.8,verbose=F)
	{
	# returns:
	#   outcome = TRUE  ==> don't reject null hypothesis; the counts are uniform
	#           = FALSE ==> reject null hypothesis;       the counts are biased
	#           = NA    ==> unable to run the test
	#   pvalue          ==> chisq.test's pvalue
	#   alpha           ==> alpha used for the test
	#   df              ==> the test's degrees of freedom
	#   sampleSize      ==> number of items

	if (any(case<5))
		{
	    return (list(outcome    = NA,
	                 pvalue     = NA,
	                 alpha      = NA,
	                 df         = NA,
	                 sampleSize = NA))
		}
	stopifnot(case>=5)
	df <- length(case)-1
	sampleSize <- sum(case)

	alpha <- 0.0000000001 # default alpha value
	tryCatch(
		# ***** need to consider whether this should be treated as "unable to run the test" *****
		# Calculating alpha value for fixed power
		alpha <- pwr.chisq.test(w=effectSize,df=df,N=sampleSize,power=power,sig.level=NULL)$sig.level
		,
		# But if an error occurs, we will stick to the default alpha value.
		# This happens when the alpha value becomes extremely small.
		error=function(error_message)
			{
			message("Uniroot problem: alpha will remain at the minimum value of 0.0000000001.")
			message(error_message)
			message("")  # creates a newline after error_message
		  	}
		)

	pvalue <- chisq.test(case)$p.value
	if (verbose)
		write(paste("alpha:",alpha,"pvalue:",pvalue,"df:",df,"sampleSize:",sampleSize), stderr())
	# (pvalue <= alpha) ==> reject null hypothesis;       the counts are biased
	# (pvalue >  alpha) ==> don't reject null hypothesis; the counts are uniform
    return (list(outcome    = pvalue>alpha,
                 pvalue     = pvalue,
                 alpha      = alpha,
                 df         = df,
                 sampleSize = sampleSize))
	}

