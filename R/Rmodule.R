
find.roots = function (f, interval, lower = min(interval), upper = max(interval), 
                       tol = .Machine$double.eps^0.2, maxiter = 1000, n = 100, ...) 
{
    x = seq(lower, upper, len = n + 1)
    y = f(x, ...)
    zeros = x[which(y == 0)]
    prods = y[1:n] * y[2:(n + 1)]
    pos = which(prods < 0)
    for (i in pos)
        zeros = c(zeros, uniroot(f, lower = x[i], upper = x[i + 1], maxiter = maxiter, tol = tol, ...)$root)
    if (length(pos) == 0)
        pos = c(1, n)
    results = list(zeros = zeros, signs = sign(y), pos = pos)
    results
}

T = function(x, a, b)
{
    tan(pi / (b - a) * (x - (a + b) / 2))
}

T.inv = function(y, a, b)
{
    (b - a) / pi * atan(y) + (a + b) / 2
}

#' Update the state vector of the correlation parameters.
#'
#' @details This function takes the current state of the chain and returns the next state. The correlation parameters are updated one at a time by way of a Metropolis-Hastings Gaussian random walk for each parameter. When the set of valid values for the proposal comprises a disconnected subset, i.e., two or more disjoint subintervals, of \eqn{(-1, 1)}{(-1, 1)}, the Apes of Wrath algorithm is used to update the parameter in question.
#'
#' @param r a \eqn{p}-vector of correlations, the current state of the Markov chain.
#' @param data an \eqn{n} x \eqn{d} matrix such that the rows are iid outcomes for the study in question.
#' @param R a \eqn{d} x \eqn{d} correlation matrix in symbolic form. The off-diagonal elements should be numbered from 2 to \eqn{p+1}.
#' @param log.f the log objective function, which must take the dataset, a correlation matrix, and perhaps additional arguments.
#' @param log.f.args additional arguments for \code{log.f}.
#' @param log.priors a list of log prior densities for the correlation parameters, each of which should accept a correlation and perhaps additional arguments.
#' @param log.priors.args a list of additional arguments for the functions in \code{log.priors}.
#' @param sigma a vector, the standard deviations of the Gaussian proposals for the \eqn{p} correlation parameters. This argument must have length 1 or length \eqn{p}. In the former case, all of the random-walk proposals have the same variance. In the latter case, the proposals have distinct variances.
#' @param n a positive integer, the number of grid points to employ in root finding. The default value is 100, but in some cases a larger value may be required to avoid missing roots of the determinant function.
#
#' @return a \eqn{p}-vector, the new state of the chain.
#'
#' @export
#'
#' @examples
#' 
#' # The following function computes HPD intervals.
#'
#' hpd = function(x, alpha = 0.05)
#' {
#'     n = length(x)
#'     m = round(n * alpha)
#'     x = sort(x)
#'     y = x[(n - m + 1):n] - x[1:m]
#'     z = min(y)
#'     k = which(y == z)[1]
#'     c(x[k], x[n - m + k])
#' }
#' 
#' # The following function computes the log likelihood.
#'
#' logL = function(data, R, args)
#' {
#'     n = nrow(data)
#'     Rinv = solve(R)
#'     detR = -0.5 * n * determinant(R, log = TRUE)$modulus
#'     qforms = -0.5 * sum(diag(data %*% Rinv %*% t(data)))
#'     f = detR + qforms
#'     if (f > 0)
#'         return(-1e6)
#'     f
#' }
#' 
#' # Use a Uniform(-1, 1) prior for each correlation.
#'
#' logP = function(r, args) dunif(r, -1, 1, log = TRUE)
#' 
#' # Build the list of priors and their arguments.
#'
#' log.priors = list(logP, logP, logP, logP, logP)
#' log.priors.args = list(0, 0, 0, 0, 0)
#' 
#' # Simulate a dataset to work with. The dataset will have 32 observations,
#' # each of length 4. The outcomes will be generated from a Gaussian copula
#' # model having t-distributed marginal distributions. Then we Gaussianize
#' # the ranks for analysis.
#'
#' n = 16
#' R = diag(1, 4, 4)
#' R[1, 2] = R[2, 1] = 2
#' R[3, 4] = R[4, 3] = 3
#' R[1, 3] = R[3, 1] = R[2, 4] = R[4, 2] = 4
#' R[1, 4] = R[4, 1] = 5
#' R[2, 3] = R[3, 2] = 6
#' r = c(-0.2, -0.2, -0.4, -0.7, 0.9)
#' block = R
#' for (j in 1:5)
#'     block[block == j + 1] = r[j]
#' blist = vector("list", n)
#' for (j in 1:n)
#'     blist[[j]] = block
#' C = t(chol(as.matrix(Matrix::bdiag(blist))))
#' set.seed(42)
#' z = as.vector(C %*% rnorm(n * 4))
#' u = pnorm(z)
#' y = qt(u, df = 3)
#' data = matrix(y, n, 4, byrow = TRUE)
#' data = matrix(qnorm(rank(data) / (n * 4 + 1)), n, 4)
#' 
#' # Simulate a sample path of length 1,000.
#'
#' m = 1000
#' r.chain = matrix(0, m, 5)
#' r.chain[1, ] = 0
#' sigma = c(1, 1, 0.25, 2, 5)   # proposal standard deviations
#' start = proc.time()
#' for (i in 2:m)
#'     r.chain[i, ] = update_R(r.chain[i - 1, ], data, R,
#'                             log.f = logL,
#'                             log.priors = log.priors,
#'                             log.priors.args = log.priors.args,
#'                             sigma = sigma,
#'                             n = 400)
#' stop = proc.time() - start
#' stop
#' stop[3] / m   # 0.001 seconds per iteration on a 3.6 GHz 10-Core Intel Core i9
#' 
#' # Now show trace plots along with the truth and the 95% HPD interval.
#'
#' dev.new()
#' plot(r.chain[, 1], type = "l")
#' abline(h = r[1], col = "orange", lwd = 3)
#' abline(h = hpd(r.chain[, 1]), col = "blue", lwd = 3)
#' 
#' dev.new()
#' plot(r.chain[, 2], type = "l")
#' abline(h = r[2], col = "orange", lwd = 3)
#' abline(h = hpd(r.chain[, 2]), col = "blue", lwd = 3)
#' 
#' dev.new()
#' plot(r.chain[, 3], type = "l")
#' abline(h = r[3], col = "orange", lwd = 3)
#' abline(h = hpd(r.chain[, 3]), col = "blue", lwd = 3)
#' 
#' dev.new()
#' plot(r.chain[, 4], type = "l")
#' abline(h = r[4], col = "orange", lwd = 3)
#' abline(h = hpd(r.chain[, 4]), col = "blue", lwd = 3)
#' 
#' dev.new()
#' plot(r.chain[, 5], type = "l")
#' abline(h = r[5], col = "orange", lwd = 3)
#' abline(h = hpd(r.chain[, 5]), col = "blue", lwd = 3)

update_R = function(r, data, R, log.f, log.f.args, log.priors, log.priors.args, sigma, n = 100)
{
    R.new = R.old = R
    m = length(r)
    if (length(sigma) == 1)
        sigma = rep(sigma, m)
    if (length(log.priors) == 1)
    {
        temp = log.priors[[1]]
        log.priors = vector("list", m)
        for (j in 1:m)
            log.priors[[j]] = temp
        temp = log.priors.args
        log.priors.args = vector("list", m)
        log.priors.args[1:m] = rep(temp, m)
    }
    for (j in 1:m)
        R.old[R == j + 1] = R.new[R == j + 1] = r[j]
    for (j in 1:m)
    {
        roots = find.roots(detR, lower = -1, upper = 1, param = j + 1, R = R, R_ = R.old, n = n)
        zeros = roots$zeros
        signs = roots$signs
        pos = roots$pos
        if (length(zeros) == 0)
        {
            int = c(-1, 1)
            eta.old = T(r[j], int[1], int[2])
            eta.new = eta.old + rnorm(1, sd = sigma[j])
            r.new = T.inv(eta.new, int[1], int[2])
        }
        else if (length(zeros) == 1)
        {
            if (signs[pos + 1] > 0)
                int = c(zeros, 1)
            else
                int = c(-1, zeros)
            eta.old = T(r[j], int[1], int[2])
            eta.new = eta.old + rnorm(1, sd = sigma[j])
            r.new = T.inv(eta.new, int[1], int[2])
        }
        else if (length(zeros) == 2)
        {
            if (signs[pos[1] + 1] > 0)
            {
                int = zeros
                eta.old = T(r[j], int[1], int[2])
                eta.new = eta.old + rnorm(1, sd = sigma[j])
                r.new = T.inv(eta.new, int[1], int[2])
            }
            else
            {
                int = c(-1, zeros, 1)
                int = matrix(int, ncol = 2, byrow = TRUE)
                r.new = apesOfWrath(r[j], int, sigma[j])
            }
        }
        else if (length(zeros) > 2)
        {
            if (signs[pos[1] + 1] > 0)
            {
                int = zeros
                if (length(zeros) %% 2 == 1)
                    int = c(int, 1)
            }
            else
            {
                int = c(-1, zeros)
                if (length(zeros) %% 2 == 0)
                    int = c(int, 1)
            }
            int = matrix(int, ncol = 2, byrow = TRUE)
            r.new = apesOfWrath(r[j], int, sigma[j])
        }
        R.new[R == j + 1] = r.new
        log.epsilon = log.f(data, R.new, log.f.args) - log.f(data, R.old, log.f.args)
        log.epsilon = log.epsilon + log(1 + eta.new^2) - log(1 + eta.old^2)
        log.epsilon = log.epsilon + log.priors[[j]](r.new, log.priors.args[[j]]) - log.priors[[j]](r[j], log.priors.args[[j]])
        if (log(runif(1)) < log.epsilon)
        {
            r[j] = r.new
            R.old[R == j + 1] = R.new[R == j + 1] = r.new
        }
        else
            R.new[R == j + 1] = r[j]
    }
    r
}
