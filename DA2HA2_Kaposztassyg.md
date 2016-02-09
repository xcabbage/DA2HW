---

title: '**DA2 Homework Assignment 2**'
output:
  word_document:
    fig_height: 4
    fig_width: 10
    
---

# NASDAQ Individual study on IBM tickets
## by Gabor Kaposztassy



```{r, echo=FALSE, warning=FALSE, message=FALSE, comment=FALSE}
library(readr)
library(dplyr)
library(dummies)
library(ggplot2)
library(lmtest)  # for Breush-Godfrey
library(stargazer)
library(pander)
library(stats)
library(data.table)


# function collection

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

pperron <- function(x, model = c('constant', 'trend'), type = "Z-tau") {
    if (!require(urca)) stop("Required urca package is missing.")

    results <- ur.pp(x, type = type, model = model)
    print(results)

    model <- match.arg(model)
    if (model == 'trend') trend = 'ct' else trend = 'c' 
    cat(
        "MacKinnon approximate p-value for Z-tau:", 
        punitroot(results@teststat, trend = trend), 
        "\n\n"
    )
}

lags <- function(variable, lags) {
    var_string <- deparse(substitute(variable))
    paste(
        lapply(
            lags, 
            function(i) {
                paste0("lag(", var_string, ",", i, ")")
            }
        ),
        collapse = "+"
    )
}

d <- function(x) {
    c(NA, diff(x))
}

Arima <- function(..., transform.pars = FALSE) {
    model <- arima(...)

    # rename to be consistent with lm
    names(model$coef) <- gsub('intercept', '(Intercept)', names(model$coef))
    row.names(model$var.coef) <- gsub('intercept', '(Intercept)', row.names(model$var.coef))
    colnames(model$var.coef) <- gsub('intercept', '(Intercept)', colnames(model$var.coef))

    model
}

capwords <- function(s, strict = FALSE) {
    cap <- function(s) paste(toupper(substring(s, 1, 1)),
                  {s <- substring(s, 2); if(strict) tolower(s) else s},
                             sep = "", collapse = "-" )
    sapply(strsplit(s, split = "-"), cap, USE.NAMES = !is.null(names(s)))
}


stargazer_r <- function(list_of_models, type="text", align=TRUE, no.space=TRUE,
                        omit.stat=c("LL", "aic", "ser", "f", "adj.rsq", "sigma2"), 
                        se = 'robust', max_lag = 0, ...) {
    if (!require(stargazer)) stop("Required stargazer package is missing.")
    
    if (class(type) != "character") stop("Different models should be given in a list.")
    if (class(list_of_models) == "lm") list_of_models <- list(list_of_models)
    if (!length(se) %in% c(1, length(list_of_models))) stop("For parameter se you should give one string (if you want to apply it to all models) or a list of strings (if you want to apply different types of standard error for the different models). The string could take traditional, robust, and newey-west (default is robust).")

    if (length(se) == 1) {
        note <- paste(capwords(se[[1]]), "standard errors in parentheses")
        se <- as.list(rep(se[[1]], length(list_of_models)))
    } else {
        note <- "Standard errors in parentheses"
    }

    if (length(max_lag) == 1) {
        max_lag <- as.list(rep(max_lag[[1]], length(list_of_models)))
        if (all(se == 'newey-west')) {
            note <- paste(note, "- max lag:", max_lag[[1]])
        }
    }

    if (any(se == 'newey-west')) omit.stat <- c(omit.stat, 'rsq')

    list_se_robust <- lapply(
        seq_along(list_of_models), 
        function(j) {
            if (class(list_of_models[[j]]) == 'lm') {
                calculate_se(list_of_models[[j]], se = se[[j]], max_lag = max_lag[[j]])
            } else {
                NULL
            }
        }
    )
    
    args <- list(...)
    if (!is.null(args[['out']])) type="html"
    
    stargazer(
        list_of_models,
        se = list_se_robust,
        report ="vcs*",
        notes = note,
        type = type, align = align, omit.stat = omit.stat, no.space = no.space,
        ...
    )
}

calculate_se <- function(lm_model, se = 'robust', max_lag) {
    if (!se %in% c('traditional', 'robust', 'newey-west')) stop("se should be one of traditional, robust or newey-west (default is robust).")
    if (!require(sandwich)) stop("Required sandwich package is missing.")

    if (se == 'robust') {
        sqrt(diag(vcovHC(lm_model, type="HC1")))
    } else if (se == 'newey-west') {
        sqrt(diag(NeweyWest(lm_model, lag = max_lag, prewhite = FALSE)))
    } else {
        sqrt(diag(vcov(lm_model)))
    }
}

summary_r <- function(model, se = 'robust', max_lag = 0, ...) {

    sumry <- summary(model)
    table <- coef(sumry)
    table[, 2] <- calculate_se(model, se, max_lag)
    table[, 3] <- table[,1]/table[, 2]
    table[, 4] <- 2*pt(abs(table[, 3]), df.residual(model), lower.tail=FALSE)

    sumry$coefficients <- table
    p <- nrow(table)
    if (p > 1) {
        if (se == 'robust') {
            hyp <- cbind(0, diag(p - 1))    
            sumry$fstatistic[1] <- linearHypothesis(model, hyp, white.adjust="hc1")[2, "F"]
        } else if (se == 'newey-west') {
            sumry$fstatistic[1] <- NA
        }
    }

    print(sumry)
    cat("Number of observations:", length(residuals(model)), "\n\n")

    if (se == 'robust') {
        cat("Note: Heteroscedasticity-consistent standard errors (adjustment HC1)\n")
    } else if (se == 'newey-west') {
        cat("Note: Newey-West standard errors - maximum lag:", max_lag, "\n")
    }
    

}

summ <- function(x) {
    cbind(mean = mean(x), median = median(x), var = var(x), sd = sd(x))
}

```

The following document contains the times series analysis of IBM daily closing prices and volumes from 1962. I used adjusted closing price due to the 7 stock cut during the whole period. IBM and NASDAQ historical data were dowloaded from: <http://finance.yahoo.com/q/hp?s=IBM+Historical+Prices>  

## I. Basic statistics and graphs

The distribution of closing price has quite normal distribution separately in the three time-interval with high variance showed below, the volume has much lower variance.

```{r, echo=FALSE}

setwd("C:/Users/KG/Documents/CEU/DA2/HW2")
nasdaq <- read_csv('NASDAQ.csv')
ibm <- read.csv('IBM.csv')
ibm$date <- as.Date(ibm$date)
ibm$volume <- as.numeric(ibm$volume)
ibm$close <- ibm$`adj.close`
ibm <- na.omit(ibm)

pander(
rbind(
  cbind.data.frame(variable='Close USD',summ(ibm$close)),
  cbind.data.frame(variable='Volume (MUSD)',summ(ibm$volume/1000000))
  )
)

```


Close price and volume over time

```{r, echo=FALSE}

g1<-ggplot(ibm, aes(x = date, y = close)) + labs(y = 'close price') + geom_line()
g2<-ggplot(ibm, aes(x = date, y = volume/1000000)) + labs(y = 'volume (mUSD)') + geom_line()
multiplot(g1,g2)

```


New variables were derived for analysis purpose:

* dummies for months and weekdays
* indentified time-gaps (911, Sandy, holidays)
* log of close price and volume and its logdifs as well
* dummies for bigger trend breaks

```{r, echo=FALSE}
# Data manipulation

Sys.setlocale("LC_TIME", "English") # Windows
ibm <- ibm %>%
    mutate(
        year = as.numeric(format(date, '%Y')),
        month = format(date, '%b'),
        day_of_week = weekdays(date)
    ) %>%
    arrange(date) %>%
    mutate(gap = factor(as.numeric(date - lag(date)) - 1))


ibm <- dummy.data.frame(as.data.frame(ibm))
names(ibm) <- gsub('month|day_of_week', '', names(ibm))

ibm <- ibm %>%
    rename(after_weekend = gap2, after_911 = gap6) %>%
    mutate(
        after_sandy = as.numeric(date == '2012-10-31'),
        after_holiday = as.numeric((gap3 > 0 | gap4 > 0) & after_sandy == 0),
        before_holiday = lead(after_holiday)
    )
ibm <- ibm %>%
    mutate(
        t1 = as.numeric(date < '1998-01-01'),
        t2 = as.numeric(date >= '1998-01-01' & date < '2008-10-01'),
        t3 = as.numeric(date >= '2008-10-01')
    )

ibm <- ibm %>%
    mutate(
        ln_close = log(close),
        ln_volume = ifelse(volume == 0, NA, log(volume)),  # to have NA instead of Inf

        lnreturn = ln_close - lag(ln_close),        
        dln_volume = ln_volume - lag(ln_volume)
    )
```


## Time series analisys

**Closing price**

At first let us inspect the closing price. It has different distribution during the three period:

```{r, echo=FALSE}

par(mfrow=c(1, 3))
hist(filter(ibm,t1 == 1)$close, main="1962-1998",xlab="Close price")
hist(filter(ibm,t2 == 1)$close, main="1998-2008",xlab="Close price")
hist(filter(ibm,t3 == 1)$close, main="2009-",xlab="Close price")

```

The unit-root test shows that it is not stationary, so I took the log and log dif. Only the log return is stationary with high confidence. But as it is like white noise - not predictable.
**`pperron(ibm$close)`**
**`pperron(ibm$ln_close)`**
**`pperron(ibm$lnreturn)`**

```{r, echo=FALSE}

# unit root tests
pperron(ibm$close)
pperron(ibm$ln_close)
pperron(ibm$lnreturn)

par(mfrow=c(1, 1))
plot(ibm$date, ibm$ln_close, type ="l", xlab = "Date", ylab = "ln close")
plot(ibm$date, ibm$lnreturn, type ="l", xlab = "Date", ylab = "log return")

```

Estimating an OLS with Newey-West SE, 2 lags on log returns shows no evidence for stock series prediction (only 911 and Sandy were significant), so the Efficient Market Hypotheses is true.

```{r, echo=FALSE}

days <- c('Monday', 'Tuesday', 'Wednesday', 'Friday')
x_vars <- c(
    'Jan', days, 
    'before_holiday', 'after_holiday', 'after_911', 'after_sandy')

myregs <- lm(
        as.formula(paste(
            'lnreturn ~', 
            lags(lnreturn, 1:2), 
            '+', paste(x_vars, collapse = '+')        
        )),
        data = ibm       
    )

stargazer_r(
    myregs, se = 'newey-west', max_lag = 2,
    dep.var.labels = 'log return'
)

```

Comparing the different periods, the result is the same, except that some inverse Monday occured between 1998-2008.

```{r, echo=FALSE}

myreg1 <- lm(
    as.formula(paste(
        'lnreturn ~', 
        lags(lnreturn, 1:2), 
        '+', paste(x_vars, collapse = '+')
    )),
    data = filter(ibm, year <= 1988)   
)

myreg2 <- lm(
    as.formula(paste(
        'lnreturn ~', 
        lags(lnreturn, 1:2), 
        '+', paste(x_vars, collapse = '+')
    )),
    data = filter(ibm, year > 1988, year <= 2008)   
)

myreg3 <- lm(
    as.formula(paste(
        'lnreturn ~', 
        lags(lnreturn, 1:2), 
        '+', paste(x_vars, collapse = '+')
    )),
    data = filter(ibm, year > 2008)   
)

stargazer_r(
    list(myreg1, myreg2, myreg3), se = 'newey-west', max_lag = 2,
    dep.var.labels = "log return"
)

```

**Volume**

Volume looks sationary based on the unit-root test, but for any case let us take log and log dif of volume:

**`pperron(ibm$volume)`**
**`pperron(ibm$ln_volume)`**
**`pperron(ibm$dln_volume)`**

```{r, echo=FALSE}
# unit root tests
pperron(ibm$volume)
pperron(ibm$ln_volume)
pperron(ibm$dln_volume)

par(mfrow=c(1, 1))
plot(ibm$date, ibm$ln_volume, type ="l", xlab = "Date", ylab = "ln volume")
plot(ibm$date, ibm$dln_volume, type ="l", xlab = "Date", ylab = "log dif of volume")

```

ACF and PACF of log dif of volume suggest 2 or 3 lags:
```{r, echo=FALSE}
# correlograms
par(mfrow=c(1, 2))
acf(na.omit(ibm$dln_volume), main='')
pacf(na.omit(ibm$dln_volume), main='')

```

Regressing the ARMA(2,2) model on log dif of volume

```{r, echo=FALSE}

arimareg <- Arima(ibm$dln_volume, order = c(2, 0, 2), xreg = ibm[x_vars])
#resid <- na.omit(arimareg$residuals)
#acf(resid)
#pacf(resid)

stargazer_r(
    list(arimareg), se = 'newey-west', max_lag = 2,
    dep.var.labels = 'ln volume'
)

```

There is no trend break in volume log dif.


## Regression with NASDAQ ##

The β of NASDAQ log return shows no significant correlation.

```{r, echo=FALSE}
nasdaq <- read_csv('NASDAQ.csv')

nasdaq <- nasdaq %>%
    mutate(
        nas_ln_close = log(close),
        nas_ln_volume = ifelse(volume == 0, NA, log(volume)),  # to have NA instead of Inf
        nas_lnreturn = nas_ln_close - lag(nas_ln_close),        
        nas_dln_volume = nas_ln_volume - lag(nas_ln_volume),
        nas_date = date
    )

nasdaq <- na.omit(nasdaq)

setDT(ibm)
setDT(nasdaq)
setkey(nasdaq, nas_date)
setkey(ibm, date)

ibm <- merge(ibm, nasdaq[,8:12, with = FALSE], by.x="date", by.y = "nas_date")

ibm <- na.omit(ibm)

md1 <- lm(lnreturn ~ nas_lnreturn, data = ibm)

x_vars <- c(
    'Jan', days, 'nas_lnreturn',
    'before_holiday', 'after_holiday', 'after_911', 'after_sandy')

myregs <- lm(
        as.formula(paste(
            'lnreturn ~', 
            lags(lnreturn, 1:2), 
            '+', paste(x_vars, collapse = '+')        
        )),
        data = ibm       
    )

stargazer_r(
    list(md1,myregs), se = 'newey-west', max_lag = 2,
    dep.var.labels = 'log return'
)

```

With ARMA(2,2) the IBM log returns shows significant correlation with the NASDAQ composite log return, 10 pp higher NASDAQ index expected a 6.9 pp higher IBM index.

```{r, echo=FALSE}
ibm <- data.frame(ibm)
arima1 <- Arima(ibm$lnreturn, order = c(2, 0, 2), xreg = ibm$nas_lnreturn)
arima2 <- Arima(ibm$lnreturn, order = c(2, 0, 2), xreg = ibm[x_vars])

stargazer_r(
    list(arima1,arima2), se = 'newey-west', max_lag = 2,
    dep.var.labels = 'log return'
)

```

