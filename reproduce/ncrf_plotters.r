# plot_ncrf_event_matrix--
#   Plot, one dot per alignment, error events vs alignment length.
#
# The input file is the output of "ncrf_extract_event_matrix --withheader".
# It has a header line and (usually) one row for each alignment;  9 columns
# per row.
#
# typical input:
#   line motif mRatio m    mm  io  ix do  dx
#   1    GGAAT 0.817  873  119 38  0  36  3
#   5    GGAAT 0.815  1665 207 68  0  92  10
#   9    GGAAT 0.816  2776 334 114 0  161 17
#   13   GGAAT 0.82   545  72  17  0  27  4
#   17   GGAAT 0.826  881  97  30  0  50  8
#    ...

plot_ncrf_event_matrix <- function(name,filename,motif=NULL,
                                   xMax=NA,yMax=NA,mainFontSize=0.9,
                                   width=5,height=5.6,newWindow=T,pdfName=NULL)
	{
	dd = read.table(filename,header=T,
	                colClasses=c("integer","character","numeric",
	                             "integer","integer","integer",
	                             "integer","integer","integer"))

	if (!is.null(motif))
		dd = dd[(dd[,"motif"]==motif),]

	dd[,"events"] = dd[,"m"]+dd[,"mm"]+dd[,"io"]+dd[,"ix"]+dd[,"do"]+dd[,"dx"]
	dd.xMax = max(dd[,"mm"],dd[,"io"]+dd[,"ix"],dd[,"do"]+dd[,"dx"])
	dd.yMax = max(dd[,"events"])

	if (is.na(xMax))
		xMax = dd.xMax
	else if (xMax < dd.xMax)
	    warning(paste("xMax=",xMax," hides some points; (actual xMax is ",dd.xMax,")",sep=""))

	if (is.na(yMax))
		yMax = dd.yMax
	else if (yMax < dd.yMax)
	    warning(paste("yMax=",yMax," hides some points; (actual yMax is ",dd.yMax,")",sep=""))

	mmSlope = sum(dd[,"mm"]) / sum(dd[,"events"])
	ioSlope = sum(dd[,"io"]) / sum(dd[,"events"])
	ixSlope = sum(dd[,"ix"]) / sum(dd[,"events"])
	doSlope = sum(dd[,"do"]) / sum(dd[,"events"])
	dxSlope = sum(dd[,"dx"]) / sum(dd[,"events"])

	title = paste("alignment events in",name)
	motifs = unique(dd[,"motif"])
	if (length(motifs) == 1)
		title = paste(title,"\nfor ",motifs[1],sep="")

	alignedBases = sum(dd[,"m"],dd[,"mm"])
	yLab = paste("total events in alignment (m+mm+i+d)",
	             "\n(",format(alignedBases,big.mark=",",trim=T)," bases aligned overall)",sep="")

	if (!is.null(pdfName))
		pdf(file=pdfName,width=width,height=height,pointsize=10)
	else if (newWindow)
		quartz(width=width,height=height)

	legText  = c(paste("mismatch ",format(100*mmSlope,digits=2),"%",sep=""),
	             paste("ins open ",format(100*ioSlope,digits=2),"%",sep=""),
	             paste("ins ext ", format(100*ixSlope,digits=2),"%",sep=""),
	             paste("del open ",format(100*doSlope,digits=2),"%",sep=""),
	             paste("del ext ", format(100*dxSlope,digits=2),"%",sep=""))
	legPch   = c(0,1,3,2,4)
	legColor = c("black","red","green","blue","purple")

	par(mar=c(4.1,5.1,3.1,1.1)) # BLTR
    par(cex.main=mainFontSize)

	plot(dd[,"mm"],dd[,"events"],pch=legPch[1],col=legColor[1],
		 xlim=c(0,xMax),ylim=c(0,yMax),
		 main=title,xlab="errors in alignment",ylab=yLab)
	points(dd[,"io"],dd[,"events"],pch=legPch[2],col=legColor[2])
	points(dd[,"ix"],dd[,"events"],pch=legPch[3],col=legColor[3])
	points(dd[,"do"],dd[,"events"],pch=legPch[4],col=legColor[4])
	points(dd[,"dx"],dd[,"events"],pch=legPch[5],col=legColor[5])
	lines(c(0,mmSlope*yMax),c(0,yMax),lty=2,col=legColor[1])
	lines(c(0,ioSlope*yMax),c(0,yMax),lty=2,col=legColor[2])
	lines(c(0,ixSlope*yMax),c(0,yMax),lty=2,col=legColor[3])
	lines(c(0,doSlope*yMax),c(0,yMax),lty=2,col=legColor[4])
	lines(c(0,dxSlope*yMax),c(0,yMax),lty=2,col=legColor[5])
	legend("bottomright",legend=legText,pch=legPch,col=legColor)

	if (!is.null(pdfName))
		dev.off()
	}

# plot_ncrf_filter_results--
#   Plot, one dot per alignment, accepted/rejected filtering results.
#
# The input file is the output of "error_nonuniformity_filter --report:matrix".
# It has NO header line and one row for each alignment;  2M+2 columns per row,
# where M is the motif length.  Column 1 is the line number of the alignment
# (in the input to error_nonuniformity_filter).  Column 2 is the outcome of the
# test.  Columns 3 thru M+2 are positional match counts for that alignment,
# and columns M+3 thru 2M+2 are the positional error counts.  We assume all
# rows have the same number of columns, i.e. that the same motif length is
# represented in all rows.
#
# typical input:
#   1  rejected     294  298  277  298  296  49  26  70  31  55
#   10 not_rejected 214  222  219  222  219  45  22  33  26  36
#   19 rejected     313  324  311  323  316  70  30  58  31  44
#   28 not_rejected 834  857  826  858  836  154 102 167 119 149
#   37 not_rejected 709  734  701  730  722  156 88  149 91  112
#    ...

plot_ncrf_filter_results <- function(name,filename,motifLength,
                                     width=5,height=5.6,newWindow=T,pdfName=NULL)
	{
	dd = read.table(filename,header=F,
		            colClasses=c("integer","character",rep("integer",2*motifLength)))

	resultCol = 2
	mCol = 3
	xCol = mCol + motifLength
	dd[,"m*"] = rowSums(dd[,mCol:(mCol+motifLength-1)])
	dd[,"x*"] = rowSums(dd[,xCol:(xCol+motifLength-1)])
	dd[,"mRatio"] = dd[,"m*"] / (dd[,"m*"]+dd[,"x*"])

	xMax = max(dd[,"m*"])
	yMax = max(dd[,"x*"])

	accepted = (dd[,resultCol]=="not_rejected")
	rejected = (dd[,resultCol]=="rejected")

	title = paste("alignment filtering in",name)
	xLab = "errors in alignment"
	yLab = "matches in alignment"

	if (!is.null(pdfName))
		pdf(file=pdfName,width=width,height=height,pointsize=10)
	else if (newWindow)
		quartz(width=width,height=height)

	legText  = c(paste("accepted (",sum(accepted),")",sep=""),
	             paste("rejected (",sum(rejected),")",sep=""))
	legPch   = c(0,1)
	legColor = c("black","red")

	dd[accepted,"pch"]   = legPch[1]
	dd[rejected,"pch"]   = legPch[2]
	dd[accepted,"color"] = legColor[1]
	dd[rejected,"color"] = legColor[2]

	shuffle = sample(nrow(dd),nrow(dd))

	par(mar=c(4.1,4.1,3.1,1.1)) # BLTR

	plot(dd[shuffle,"m*"],dd[shuffle,"x*"],
	     pch=dd[shuffle,"pch"],col=dd[shuffle,"color"],
	     xlim=c(0,xMax),ylim=c(0,yMax),
		 main=title,xlab=xLab,ylab=yLab)
	legend("bottomright",legend=legText,pch=legPch,col=legColor)

	if (!is.null(pdfName))
		dev.off()
	}

# read_ncrf_summary--
#   Read the output of ncrf_summary.py.
#
# typical input:
#   #line  motif seq               start end   strand seqLen querybp mRatio m     mm   i   d
#   1      GGAAT SRR2036394.100749 9566  10445 +      14279  823     81.7%  730   78   71  15
#   10     GGAAT SRR2036394.100749 11581 12099 +      14279  505     83.8%  446   45   27  14
#   19     GGAAT SRR2036394.101301 12191 12923 +      47491  744     83.0%  635   76   21  33
#    ...

read_ncrf_summary <- function(filename,motif=NULL)
	{
	dd = read.table(filename,header=T,comment.ch="",
	                colClasses=c("integer","character","character",
	                             "integer","integer","character",
	                             "integer","integer","character",
	                             "integer","integer","integer","integer"))
    colnames(dd)[1] = "line"

	if (!is.null(motif))
		dd = dd[(dd[,"motif"]==motif),]

    dd$mRatio <- as.numeric(gsub("%","",dd$mRatio))/100

	dd
	}
