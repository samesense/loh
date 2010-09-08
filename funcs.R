library(ggplot2)
library(cnv)

plot.cnv.all.perry <- function(data, chrom.gap=2e7, colour=5, title=NA, ylim=c(-2,2), xlabel='Chromosome')
{
	# modified from cnv-seq

	#if(nrow(subset(data, data$log2>max(ylim)|data$log2<min(ylim)))>0) warning('missed some data points due to small ylim range')
	level <- sort(as.character(unique(data$chromosome)))
	ok.to.convert <- !is.na(suppressWarnings(as.numeric(level)))
	level <- c(sort(as.numeric(level[ok.to.convert])), level[!ok.to.convert])
	labels<-breaks<-c()
	last.pos <- 0
	temp <- data.frame()
	for(chr in level)
	{
		sub <- subset(data, chromosome==chr)
		if(nrow(sub)==0) next
		sub$gpos <- last.pos+sub$position
		labels <- c(labels, chr)
		breaks <- c(breaks, (min(sub$gpos)+max(sub$gpos))/2)
		last.pos <- max(sub$gpos)+chrom.gap
		temp <- rbind(temp, sub)
	}
	data <- temp
	data$colour <- as.character(match(data$chromosome,level)%%colour)
	p <- ggplot()+geom_point(data=data, map=aes(x=gpos, y=log2, colour=colour))+opts(legend.position='none')
	p <- p + scale_y_continuous(expression(paste(Log[2],' Ratio')))
	p <- p + scale_x_continuous(xlabel, breaks=breaks, labels=labels)
#	cat(breaks,"\n", sep=',')
#	cat(labels,"\n",sep=',')
	if(colour==2) p <- p + scale_colour_discrete(l=30)
	else p <- p + scale_colour_brewer(palette='Set1')
	p$legend.position <- 'none'
	if(!is.na(title)) p$title <- title
	p
}

plot.murim <- function(data, chrom.gap=2e7, colour=5, title=NA, ylim=c(-2,2), xlabel='Chromosome')
{
	# Murim's plot
	level <- sort(as.character(unique(data$CHR)))
	ok.to.convert <- !is.na(suppressWarnings(as.numeric(level)))
	level <- c(sort(as.numeric(level[ok.to.convert])), level[!ok.to.convert])
	labels<-breaks<-c()
	last.pos <- 0
	temp <- data.frame()
	for(chrsme in level)
	{
		sub <- subset(data, CHR==chrsme)
		if(nrow(sub)==0) next
		sub$gpos <- last.pos+sub$MapInfo
		labels <- c(labels, chrsme)
		breaks <- c(breaks, (min(sub$gpos)+max(sub$gpos))/2)
		last.pos <- max(sub$gpos)+chrom.gap
		temp <- rbind(temp, sub)
	}
	data <- temp
	data$colour <- as.character(match(data$CHR,level)%%colour)
	p <- ggplot()+geom_point(data=data, map=aes(x=gpos, y=Diff, colour=colour))+opts(legend.position='none')
	p <- p + scale_y_continuous(expression('Allele Diff'))
	p <- p + scale_x_continuous(xlabel, breaks=breaks, labels=labels)
#	cat(breaks,"\n", sep=',')
#	cat(labels,"\n",sep=',')
	if(colour==2) p <- p + scale_colour_discrete(l=30)
	else p <- p + scale_colour_brewer(palette='Set1')
	p$legend.position <- 'none'
	if(!is.na(title)) p$title <- title
	p
}