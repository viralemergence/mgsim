griwm <- function(date, sd.vec, iter, extinv, verbose = FALSE) {
  # Transform if first arrival instead of extinction
  k <- length(date)
  iter <- 1000; itdiv <- iter/(iter/10)
  T.up.vec <- T.mci.vec <- w.T.mci.vec <- rep(0,iter)
  T.up.vec <- T.mci.vec <- w.T.mci.vec <- rep(0,iter)

  for (c in 1:iter) {
    date.samp <- rep(0,k)

    for (b in 1:k) {
      date.samp[b] <- round(rnorm(1,date[b],sd.vec[b]))}

    if (extinv == "extinction") {
      date.samp <- (sort(date.samp))
      date.samp2 <- max(date.samp) + 1 - date.samp}

    if (extinv == "invasion") {
      date.samp <- rev(sort(date.samp))
      date.samp2 <- min(date.samp) - 1 + date.samp}

    k <- length(date.samp2)
    e <- matrix(1,nrow=k,ncol=1)
    dat.rnge <- date.samp2[1] - date.samp2[k]
    dat.diff <- date.samp2[1] - date.samp2[2:k]
    dat.diff <- dat.diff[dat.diff > 0] ## removes zeros if two duplicate end dates
    nu <- (1/(k-1)) * sum(log(dat.rnge/na.omit(dat.diff[1:(k-2)])))

    ## k x k matrix
    kmat <- matrix(0, nrow=k, ncol=k)

    for (j in 1:k) {
      for (i in 1:k) {
        kmat[i,j] <- ifelse(j <= i, (gamma((2*nu) + i) * gamma(nu+j)) / (gamma(nu+1) * gamma(j)), 0)}}   ## for j <= i

    llq.kmat <- kmat
    urq.kmat <- t(kmat)
    final.kmat <- llq.kmat + urq.kmat   # creates symmetric kmat
    diag(final.kmat) <- diag(kmat)

    w <- MASS::ginv(t(e) %*% MASS::ginv(final.kmat) %*% e) %x% (MASS::ginv(final.kmat) %*% e)

    if (extinv == "extinction") {
      theta <- max(date.samp) + 1 - sum(w*date.samp2)}
    if (extinv == "invasion") {
      theta <- -min(date.samp) + 1 + sum(w*date.samp2)}

    alpha <- 0.5
    S.lo <- (-log((1 - alpha/2))/k)^-nu
    S.up <- ((-log(alpha/2))/k)^-nu

    if (extinv == "extinction") {
      T.lo <- max(date.samp) + 1 - (date.samp2[1] + ((date.samp2[1] - date.samp2[k])/(S.lo - 1)))
      T.up <- max(date.samp) + 1 - (date.samp2[1] + ((date.samp2[1] - date.samp2[k])/(S.up - 1)))}
    if (extinv == "invasion") {
      T.lo <- -min(date.samp) + 1 + (date.samp2[1] + ((date.samp2[1] - date.samp2[k])/(S.lo - 1)))
      T.up <- -min(date.samp) + 1 + (date.samp2[1] + ((date.samp2[1] - date.samp2[k])/(S.up - 1)))}
    T.up.vec[c] <- round(T.up,0)

    ## apply McInerny et al 2006 Conserv Biol algorithm
    if (extinv == "extinction") {
      date.mci <- rev(date.samp2)}
    if (extinv == "invasion") {
      date.samp.b <- -date.samp + max(date.samp) + 1
      date.mci <- rev(max(date.samp.b) + 1 - date.samp.b)}

    t.n <- date.mci[k]
    n <- k

    i <- t.n; p.iter <- 1
    while(p.iter > alpha)
    {
      i <- i + 1
      p.iter <- (1 - (n/t.n))^(i - t.n)
    }

    if (extinv == "extinction") {
      T.mci <- max(date.samp) + 1 - i}
    if (extinv == "invasion") {
      T.mci <- max(date.samp.b) + 1 - i
      T.mci <- max(date.samp) + 1 - T.mci}

    T.mci.vec[c] <- round(T.mci,0)

    ## calculate weighted McInerney date & confidence interval
    last.diff <- 1/(date.samp-date.samp[1])[-1]
    weight <- last.diff/last.diff[1]

    if (last.diff[1] == Inf) {
      weight <- last.diff/last.diff[2]
      weight <- weight[-1]}

    ldate <- length(date.samp)
    T.mci.lst.vec <- rep(0,ldate-1)

    for (m in 1:(ldate-1)) {
      date.it <- date.samp[1:(1+m)]
      date.age.it <- date.samp[1:(1+m)]

      if (extinv == "extinction") {
        date.mci.it <- rev(max(date.it) + 1 - date.it)}
      if (extinv == "invasion") {
        date.samp.b <- -date.samp[1:(1+m)] + max(date.samp[1:(1+m)]) + 1
        date.mci.it <- rev(max(date.samp.b) + 1 - date.samp.b)}

      k <- length(date.it)
      t.n <- date.mci.it[k]
      n <- k
      T.rng <- t.n - date.mci.it[1]

      i <- t.n; p.iter <- 1
      while(p.iter > alpha)
      {
        i <- i + 1
        p.iter <- (1 - (n/t.n))^(i - t.n)
      }

      if (extinv == "extinction") {
        T.mci.lst.vec[m] <- max(date.it) + 1 - i}

      if (extinv == "invasion") {
        T.mci.lst.vec[m] <- max(date.samp.b) + 1 - i
        T.mci.lst.vec[m] <- max(date.samp) + 1 - T.mci.lst.vec[m]}
    }

    if (last.diff[1] == Inf) {
      w.T.mci.vec[c] <- round((sum(weight*T.mci.lst.vec[-1]))/sum(weight),0)}

    if (last.diff[1] != Inf) {
      w.T.mci.vec[c] <- round((sum(weight*T.mci.lst.vec))/sum(weight),0)}

    w.T.mci.vec[c]
    if (verbose) {
      if(c%%itdiv==0) print(c)
    }
  }

  prb <- 0.5
  T.up.vec.lo <- quantile(T.up.vec,probs=(1-prb/2)); T.up.vec.med <- median(T.up.vec,na.rm=T); T.up.vec.up <- quantile(T.up.vec,probs=(prb/2))
  T.mci.vec.lo <- quantile(T.mci.vec,probs=(1-prb/2)); T.mci.vec.med <- median(T.mci.vec,na.rm=T); T.mci.vec.up <- quantile(T.mci.vec,probs=(prb/2))
  T.wmci.vec.lo <- quantile(na.omit(w.T.mci.vec),probs=(1-prb/2)); T.wmci.vec.med <- median(na.omit(w.T.mci.vec)); T.wmci.vec.up <- quantile(na.omit(w.T.mci.vec),probs=(prb/2))

  ## Summary
  return(c(low = round(T.wmci.vec.lo, 0), median = round(T.wmci.vec.med, 0), high = round(T.wmci.vec.up, 0)))
}
