\name{SmoothWin}
\alias{SmoothWin}
\title{
Implementation of the soft windowing for linear models
}
\description{
  Implementation of the soft windowing on a range of methods by imposing weights on the model. The parameters k and l control the shape and bandwidth of the windowing function. 
}
\usage{
SmoothWin(object                                          ,
         data                                             ,
         t                                                ,
         m                                                ,
         l = function(ignore.me.in.default) {
           r = SmoothWin:::lseq(
             from = 1                                          ,
             to = max(abs(t[m] - min(t)), abs(t[m] - max(t)),1),
             length.out = min(500, max(1,diff(range(t))))
           )
           r = unique(round(r))
           return(r)
         }                                                ,
         k = SmoothWin:::lseq(from = .5                   ,
                              to = 10                     ,
                              length.out = 50)            ,
         min.obs   = function(ignore.me.in.default) {
           lutm = length(unique(t[m]))
           r = ifelse(lutm > 1, 35, max(pi * sqrt(length(t)), 35))
           r = max(r * lutm, length(m), na.rm = TRUE)
           r = min(r, length(t), na.rm = TRUE)
           return(r)
         }                                                ,
         weightFUN = function(x) {
           x
         }                                                ,
         residFun = function(x) {
                       resid(x)
                     }                                    ,
         predictFun = function(x) {
                       predict(x)
                   }                                      ,
         weightORthreshold = 'weight'                     ,
         cdf = plogis                                     ,
         check = 2                                        ,
         sensitivity = c(1, 1, 1)                         ,
         pvalThreshold = c(0,0,0)                         ,
         threshold = sqrt(.Machine$double.eps)*10         ,
         zeroCompensation = 0                             ,
         messages  = TRUE                                 ,
         seed      = NULL                                 ,
         simple.output = FALSE                            ,
         debug         = FALSE                            ,
         ...)
    }
\arguments{
  \item{object}{
  The fitted model. The object must support "update(weights =)". See examples
}
  \item{data}{
  data.frame. Input data that is used to fit the initial model
}
  \item{t}{
  Vector of (continues) time values.
}
  \item{m}{
  Vector of integers (peaks). Mode indices on the time component. For example 10, 11, 12. Note that it is different from t[10], t[11], t[12]
}
  \item{l}{
  Vector of continues values. Bandwidth parameter. The default uses the maximum distance of the modes (t[m]) from the time boundaries, max(max(t)-t[m],t[m]-min(t)) split on 500 points on the logarithmic scale.
}
  \item{k}{
  Vector of continues values. Shape parameter. The default uses 50 split of the values from 0.5 to 10 on the logarithmic scale.
}
  \item{min.obs}{
  Single value. The minimum observations that must be in the window. The default uses the following steps.\cr
    1. If there are more than one modes in the data, then: 35*(the number of modes)\cr
    2. If there is a single mode in the data, then: max(pi * sqrt(length(t)), 35)\cr
  * min.obs must not be less than the number of observations on the mode times\cr
  ** min.obs must be less than the total observations
}
  \item{weightFUN}{
  Weight function. By default, a vector of weights call "ModelWeight" is passed to this function. See the examples. 
}
  \item{residFun}{
  Residual computation function. The default is resid(). However, the the user can define its own function. Note that the input of the function is the model object. The default is residFun = function(object)\{resid(object)\}
  }
  \item{predictFun}{
   Similar to residFun but instead defines the predict() function. The default is predictFun = function(object)\{predict(object)\}
  }
  \item{weightORthreshold}{
  select between 'weight' (default) or 'threshold'. If set to 'weight' then the sum of weights (active sample score (ASS)) would be used as the total number of (active) observations in the window, otherwise, total weights that are greater than the threshold (count weights>=threshold) would be used for the total number of samples in the window (see `threshold`).
  }
  \item{cdf}{
  A cdf function preferably symmetric. The cdf function is used for the (window) weight generating function. The function must have two parameters precisely a location such as mean and a scale. Standard cdf functions such as pnorm, pcauchy and plogis (default) can be used. For an example of custom made function we define uniform function as below:\cr\cr
  punif0 = function(x,mean=.5,sd=sqrt(1/12))\{\cr
         a = mean - sqrt(3) *sd\cr
         b = mean + sqrt(3) *sd\cr
      r =  punif(q = x,min = a,max = b) \cr
      return(r)\cr
    \}
  }
  \item{check}{
  Single integer in \{0,1,2\}. \cr
  check=1, the function selects the times (t) with more than one observations. Further, the function only selects the values with weights greater than the threshold parameter. Mostly useful in fitting linear mixed model\cr
  check=2 (default), the function only selects the values with weights greater than the `threshold` parameter\cr 
  check=0, disables all checks
}
  \item{sensitivity}{
  Vector of three values (m,v,m*v). For example (default) c(1,1,1) specifies the same weights for mean, variance and mean*variance interaction in determining the optimal (final) window. We should stress that the window size is calculated by detecting the changes in mean and variance of predictFun() and residFun() respectively.
  }  
  \item{pvalThreshold}{
  Vector of three values. It would be used as the significant level for the mean and variation tests. If all zero (default) or negative (<= 0) then the internal adaptive method (sensitivity - see above) would be used.
  }
  \item{threshold}{
  Single positive value. The minimum value for weights before removing the corresponded samples, given check=1 or check=2 and also in weightORthreshold. Default sqrt(.Machine$double.eps)*10 ~ 10^-7
}
  \item{zeroCompensation}{
   Single non-negative value. Setting to any non-negative value would replace all (weights =< zeroCompensation) with zeroCompensation. Useful for algorithms that have difficulties with zero. Default 0.
  }
  \item{messages}{
  Logical value. Set to TRUE (default) to see the messages and errors 
}
  \item{seed}{
  seed. Default NULL
}
  \item{simple.output}{
  Logical flag. Setting to TRUE leads to not exporting the list of models for l and k. Useful for preventing memory overflow. Default FALSE
}
  \item{debug}{
  Logical flag. Setting to TRUE will show some plots for the parameter selection step. Useful for debuging. Default FALSE
  }
  \item{\dots}{
  Other parameters that can be passed to the weightFUN()
}
}


\value{
  \item{final.k, final.l }{Final values for k and l}
  \item{finalModel, model.l, model.k }{List of models for l, k and the final model. }
  \item{others}{The input parameters such as x, y, t and so on}
}

\author{
Hamed Haselimashhadi <hamedhm@ebi.ac.uk>
}


\seealso{
  \code{\link{expWeight}}
}

\examples{
#################### Example in the manuscript
set.seed(1234)
par(mfrow = c(3, 1))
# Simulating data
n = 60
t = 1:n
sd = 1
m = n / 2
x = t
y = c(0 * x[t <= n / 3]                        ,
      x[t < 2 * n / 3 & t > n / 3] * 1         ,
      0 * x[t >= 2 * n / 3]) + rnorm(n, 0, sd)
# True weights
w = weights = expWeight(
  t = t                                        ,
  k = 5                                        ,
  l = n/6                                      ,
  m = m                                        ,
  plot = 0
)
# Fitting and ploting data and models
l = lm(y ~ x, weights = w)
plot(
  x                                             ,
  y                                             ,
  ylim = c(min(y), max(y) * 1.5)                ,
  col = t \%in\% seq(n / 3 + 1, 2 * n / 3 - 1) + 1,
  cex = 1.5                                     ,
  pch = 16                                      ,
  xlab = 'Time'                                 ,
  main = 'Simulated data'
)
abline(v = x[c(n / 3 + 1, 2 * n / 3 - 1)],
       lty = 2    ,
       lwd = 4    ,
       col = 'gray')
abline(l, col = 2 , lty = 2, lwd = 4)
abline(lm(y ~ x)  ,
       col = 3    ,
       lty = 3    ,
       lwd = 4)
plot(
  t,
  w ,
  type = 'b'         ,
  main = 'True weights',
  ylab = 'Weight'    ,
  xlab = 'Time'
)
# Fitting the Windowing model
r = SmoothWin(
  object = l                     ,
  data = data.frame(y = y, x = x),
  t = t                          ,
  m = m                          ,
  min.obs = 4                    ,
  debug = FALSE
)
# Plot fitted (windowed) model
plot(r, main = 'Estimated weights from WGF')


#################### Other examples
# All examples use the Orthodont dataset in the nlme package
library(nlme)
# Sort the data on the time component (age)
Orthodont =  Orthodont[order(Orthodont$age), ]
# Modes
mode = which(Orthodont$age  \%in\%  c(12))
# Time component
time = Orthodont$age
f  = formula(distance ~ Sex)


#################### Examples ####################
### Example 1. Linear model
# Method 1 (recommanded)
lm = do.call('lm', list(formula = f, data = Orthodont))
rm(f)

# Method 2 (can cause error if you pass the formula to the lm function)
# lm = lm(distance ~ Sex, data = Orthodont)

lm.result = SmoothWin(
  object = lm,
  data = Orthodont,
  t = time,
  m = mode,
  check = 0,
  weightFUN = function(x) {
    x
  },
  debug = TRUE
)
plot(
  lm.result,
  col = Orthodont$Sex,
  pch = as.integer(Orthodont$Sex),
  main = 'Simple liner model'
)


#### Example 2. Linear Model Using Generalized Least Squares
# Method 1 (recommanded)
f   = formula(distance ~ Sex)
gls = do.call('gls', list(model = f, data = Orthodont))
rm(f)

# Method 2 (can cause error if you pass the formula to the gls function)
# gls = gls(distance ~ Sex, data = Orthodont)

gls.result = SmoothWin(
  object = gls,
  data = Orthodont,
  t = time,
  m = mode,
  check = 2,
  weightFUN = function(ignore.me) {
    varFixed(~ 1 / ModelWeight) #nlme package uses the inverse weights
  },
  debug = TRUE
)
plot(
  gls.result,
  col = Orthodont$Sex,
  pch = as.integer(Orthodont$Sex),
  main = 'Linear model using GLS'
)

#### Example 3. Linear mixed model
# Method 1 (recommanded)
fixed  = formula(distance ~ Sex)
random = formula(~ 1 | Subject)
lme = do.call('lme', list(
  fixed = fixed,
  random = random,
  data = Orthodont
))
rm(fixed, random)

# Method 2 (can cause error if you pass the formula to the lme function)
# lme = lme(fixed = distance ~ Sex, random=~1|Subject , data = Orthodont)

lme.result = SmoothWin(
  object = lme,
  data = Orthodont,
  t = time,
  m = mode,
  # Remove zero weights as well as single observation dates
  check = 1,
  weightFUN = function(ignore.me) {
    varFixed(~ 1 / ModelWeight)
  },
  debug = TRUE
)
plot(
  lme.result,
  col = Orthodont$Sex,
  pch = as.integer(Orthodont$Sex),
  main = 'Linear mixed model'
)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{"Windowing", "SmoothWin" }% use one of  RShowDoc("KEYWORDS")
