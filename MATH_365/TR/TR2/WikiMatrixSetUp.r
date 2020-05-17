# NOTE: you will need to get the right path name to your copy of this file
W = read.csv("C:/Users/lziegel1/Google Drive/Teaching/Math Comp 365/Homework/TechnicalReports/TR2/WikiLinks.csv")
print(head(W))
from=W$from
to=W$to
sites = sort(union(unique(to),unique(from)))
n=length(sites)

# The following functions convert between the site labels and the matrix index
SiteToIndex = function(s) {match(s,sites)}
IndexToSite = function(s) {sites[s]}

ii = rep(0,length(to))
jj = rep(0,length(from))
for (i in 1:length(from)) ii[i] = SiteToIndex(from[i])
for (i in 1:length(to)) jj[i] = SiteToIndex(to[i])
xx = rep(1,length(jj))
library("Matrix")
A = spMatrix(nrow=n,ncol=n,i=ii,j=jj,x=xx)