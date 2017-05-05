\name{norm}
\alias{norm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Data normalization
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
norm(data, normalization = "MEAN")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{normalization}{
%%     ~~Describe \code{normalization} here~~
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
function (data, normalization = "MEAN") 
{
    if (!is.data.frame(data) & !is.matrix(data)) {
        print("the data must be a data.frame or a matrix")
    }
    normalization = toupper(normalization)
    if (identical(normalization, "MEAN")) {
        data <- scale(log2(data))
    }
    else if (identical(normalization, "QUANTILE")) {
        data <- log2(data)
        nomCol <- colnames(data)
        data <- as.data.frame(normalize.quantiles(as.matrix(data)), 
            row.names = row.names(data))
        colnames(data) <- nomCol
    }
    else if (identical(normalization, "MEDIAN")) {
        data <- log2(data)
        for (x in 1:length(data)) {
            data[, x] <- (data[, x] - median(data[, x], na.rm = T))/(mad(data[, 
                x], na.rm = T))
        }
    }
    else {
        stop("normalization must be : 'MEAN', 'QUANTILE' or 'MEDIAN'")
    }
    return(data)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
