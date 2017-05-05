\name{volcanoPlot}
\alias{volcanoPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Representing the data with volcano plots
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
volcanoPlot(data, save = FALSE, fc = 0.5, p = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{save}{
%%     ~~Describe \code{save} here~~
}
  \item{fc}{
%%     ~~Describe \code{fc} here~~
}
  \item{p}{
%%     ~~Describe \code{p} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data, save = FALSE, fc = 0.5, p = 0.05) 
{
    require(ggplot2)
    require(ggrepel)
    require(svglite)
    for (i in 1:length(data)) {
        plotTitle <- substr(colnames(data[[i]])[2], 9, nchar(colnames(data[[i]])[2]))
        values <- as.data.frame(data[[i]])
        forplot <- data.frame(x = as.numeric(values[, 4]), y = -log10(values[, 
            3]), id = as.character(values[, 1]))
        tmp <- forplot[as.numeric(forplot$y) >= -log10(fc) & 
            abs(forplot$x) > fc, ]
        p <- ggplot(forplot) + geom_point(aes(x, y, label = id, 
            color = ifelse(y >= -log10(fc) & abs(x) >= fc, "not signi", 
                "FC")), show.legend = F) + scale_color_manual(values = c("blue", 
            "red")) + geom_text_repel(data = subset(forplot, 
            abs(forplot$x) >= fc & forplot$y >= -log10(fc)), 
            aes(x, y, label = id), size = 2) + geom_vline(xintercept = fc) + 
            geom_vline(xintercept = -fc) + geom_hline(yintercept = -log10(fc)) + 
            labs(title = plotTitle, x = "log2(Fold-change)", 
                y = "-log10(P.Value)") + theme_bw()
        print(p)
        if (save) {
            ggsave(paste0("VolcanoPlot-", plotTitle, ".svg"), 
                plot = p)
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
