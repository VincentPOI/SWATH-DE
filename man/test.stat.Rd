\name{test.stat}
\alias{test.stat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
test.stat(data, stat = "t.test", exp.design, contrast)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{stat}{
%%     ~~Describe \code{stat} here~~
}
  \item{exp.design}{
%%     ~~Describe \code{exp.design} here~~
}
  \item{contrast}{
%%     ~~Describe \code{contrast} here~~
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
function (data, stat = "t.test", exp.design, contrast) 
{
    if (!is.data.frame(data) & !is.matrix(data)) {
        print("the data must be a data.frame or a matrix")
    }
    if (!is.data.frame(exp.design) & !is.matrix(exp.design)) {
        print("the experience design must be a matrix generated from model.matrix")
    }
    if (!is.data.frame(contrast) & !is.matrix(contrast)) {
        print("the contrast matrix must be  a matrix generated from makeContrasts()")
    }
    nbComp <- ncol(contrast)
    fc <- list()
    p.value <- list()
    adjust.p.value <- list()
    listResults <- list()
    for (i in 1:nbComp) {
        op <- row.names(contrast)[contrast[, i] == 1]
        ref <- row.names(contrast)[contrast[, i] == -1]
        samplesop <- row.names(design)[row(as.matrix(design[, 
            op]))[design[, op] == 1]]
        samplesref <- row.names(design)[row(as.matrix(design[, 
            ref]))[design[, ref] == 1]]
        colop <- which(colnames(data) \%in\% samplesop)
        colref <- which(colnames(data) \%in\% samplesref)
        p.value[[i]] <- apply(data, 1, function(x) {
            t.test(as.numeric(x[colref]), as.numeric(x[colop]), 
                alternative = "t")$p.value
        })
        adjust.p.value[[i]] <- p.adjust(p.value[[i]], method = "BH")
        fc[[i]] <- rowMeans(data[, colop]) - rowMeans(data[, 
            colref])
        listResults[[i]] <- data.frame(protein = c(row.names(data)))
        listResults[[i]][paste("p.value.", colnames(contrast)[i], 
            sep = "")] = c(p.value[[i]])
        listResults[[i]][paste("adjust.p.value.", colnames(contrast)[i], 
            sep = "")] = c(adjust.p.value[[i]])
        listResults[[i]][paste("fc.", colnames(contrast)[i], 
            sep = "")] = c(fc[[i]])
    }
    return(listResults)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
