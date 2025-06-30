' Correction for pleiotropy
#' 
#' This function penalyze the regulatory interactions based on pleiotropy analysis
#' 
#' @param ss Named vector containing the gene expression signature
#' @param nes Named vector containing the normalized enrichment scores
#' @param regul Regulon object
#' @param regulators Number indicating the number of top regulators to consider for the analysis or the p-value threshold for considering significant regulators
#' @param shadow Number indicating the p-value threshold for considering a significant shadow effect
#' @param targets Integer indicating the minimal number of overlaping targets to consider a pair of regulators for pleiotropy analysis
#' @param penalty Number higher than 1 indicating the penalty for the pleiotropic interactions. 1 = no penalty
#' @param method Character string indicating the method to use for computing the pleiotropy, either absolute or adaptive
#' @return Corrected regulon object

shadowRegulon <- function(ss, nes, regul, regulators=.05, shadow=.05, targets=10, penalty=2, method=c("absolute", "adaptive")) {
    method <- match.arg(method)
    pval <- pnorm(abs(nes), lower.tail=FALSE)*2
    if (regulators<1) tfs <- names(pval)[pval<regulators]
    else tfs <- names(pval)[order(pval)[1:regulators]]
    pos <- grep("--", tfs)
    if (length(pos)>0) tfs <- tfs[-pos]
    if (length(tfs)<2) return(NULL)
    tmp <- lapply(unique(tfs), function(tf1, tfs, regul, ss, nes, targets) {
        reg <- lapply(tfs[tfs != tf1], function(tf2, regul, tf1) {
            pos <- names(regul[[tf1]]$tfmode) %in% names(regul[[tf2]]$tfmode)
            list(tfmode=regul[[tf1]]$tfmode[pos], likelihood=regul[[tf1]]$likelihood[pos])
        }, regul=regul, tf1=tf1)
        names(reg) <- tfs[tfs != tf1]
        pos <- which(names(ss) %in% names(regul[[tf1]]$tfmode))
        s2 <- rank(ss[pos])/(length(ss[pos])+1)*2-1
        s1 <- abs(s2)*2-1
        s1 <- s1+(1-max(s1))/2
        s1 <- qnorm(s1/2+.5)
        tmp <- sign(nes[tf1])
        if (tmp==0) tmp <- 1
        s2 <- qnorm(s2/2+.5)*tmp
        tmp <- sapply(reg, function(x, s1, s2, targets) {
            if (length(x$tfmode)<targets) return(NA)
            pos <- match(names(x$tfmode), names(s1))
            sum1 <- sum(x$tfmode * x$likelihood * s2[pos])
            ss <- sign(sum1)
            ss[ss==0] <- 1
            sum2 <- sum((1-abs(x$tfmode)) * x$likelihood * s1[pos])
            ww <- x$likelihood/max(x$likelihood)
            return((abs(sum1) + sum2*(sum2>0)) / sum(x$likelihood) * sign(ss) * sqrt(sum(ww^2)))
        }, s1=s1, s2=s2, targets=targets)
        return(pnorm(tmp, lower.tail=FALSE))
    }, tfs=tfs, regul=regul, ss=ss, nes=nes, targets=targets)
    names(tmp) <- unique(tfs)
    pval <- unlist(tmp, use.names=F)
    names(pval) <- paste(rep(names(tmp), sapply(tmp, length)), unlist(lapply(tmp, names), use.names=FALSE), sep=" x ")
    pval <- pval[!is.na(pval)]
    regind <- t(combn(tfs, 2))
    regind <- filterRowMatrix(regind, paste(regind[, 1], regind[, 2], sep=" x ") %in% names(pval))
    pval <- cbind(pval[match(paste(regind[, 1], regind[, 2], sep=" x "), names(pval))], pval[match(paste(regind[, 2], regind[, 1], sep=" x "), names(pval))])
    tests <- table(as.vector(regind))
    switch(method,
    absolute={
        tmp <- rbind(regind[pval[, 1]<shadow & pval[, 2]>shadow, ], regind[, 2:1][pval[, 1]>shadow & pval[, 2]<shadow, ])
        if (nrow(tmp)==0) return(NULL)
        for (i in 1:nrow(tmp)) {
            ll <- regul[[tmp[i, 1]]]$likelihood
            pos <- which(names(regul[[tmp[i, 1]]]$tfmode) %in% names(regul[[tmp[i, 2]]]$tfmode))
            ll[pos] <- ll[pos]/penalty^(1/tests[tmp[i, 1]])
            regul[[tmp[i, 1]]]$likelihood <- ll
        }
    },
    adaptive={
        pval1 <- log10(pval[, 2])-log10(pval[, 1])
        tmp <- NULL
        if (length(which(pval1>0))>0) tmp <- filterRowMatrix(regind, pval1>0)
        if (length(which(pval1<0))>0) tmp <- rbind(tmp, filterRowMatrix(regind, pval1<0)[, 2:1])
        pval1 <- c(pval1[pval1>0], -pval1[pval1<0])
        if (is.null(nrow(tmp))) return(NULL)
        for (i in 1:nrow(tmp)) {
            ll <- regul[[tmp[i, 1]]]$likelihood
            pos <- which(names(regul[[tmp[i, 1]]]$tfmode) %in% names(regul[[tmp[i, 2]]]$tfmode))
            ll[pos] <- ll[pos]/(1+pval1[i])^(penalty/tests[tmp[i, 1]])
            regul[[tmp[i, 1]]]$likelihood <- ll
        }
    })
    return(regul[which(names(regul) %in% tmp[, 1])])
}