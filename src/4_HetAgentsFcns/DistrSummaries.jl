
@doc raw"""
    distrSummaries(distr,c_a_star,c_n_star,n_par,inc,incgross,m_par)

Compute distributional summary statistics, e.g. Gini indexes, top-10%
income and wealth shares, and 10%, 50%, and 90%-consumption quantiles.

# Arguments
- `distr`: joint distribution over bonds, capital and income ``(m \times k \times y)``
- `c_a_star`,`c_n_star`: optimal consumption policies with [`a`] or without [`n`]
    capital adjustment
- `n_par::NumericalParameters`, `m_par::ModelParameters`
- `inc`: vector of (on grid-)incomes, consisting of labor income (scaled by ``\frac{\gamma-\tau^P}{1+\gamma}``, plus labor union-profits),
    rental income, liquid asset income, capital liquidation income,
    labor income (scaled by ``\frac{1-\tau^P}{1+\gamma}``, without labor union-profits),
    and labor income (without scaling or labor union-profits)
- `incgross`: vector of (on grid-) *pre-tax* incomes, consisting of
    labor income (without scaling, plus labor union-profits), rental income,
    liquid asset income, capital liquidation income,
    labor income (without scaling or labor union-profits)
"""
function distrSummaries(
    distr::AbstractArray,
    Q::Real,
    c_a_star::AbstractArray,
    c_n_star::AbstractArray,
    n_par::NumericalParameters,
    inc::AbstractArray,
    incgross::AbstractArray,
    m_par::ModelParameters,
)
    ## Distributional summaries
    distr_m = sum(distr, dims = (2, 3))[:]
    distr_k = sum(distr, dims = (1, 3))[:]
    distr_y = sum(distr, dims = (1, 2))[:]

    total_wealth = Array{eltype(distr)}(undef, n_par.nk .* n_par.nm)
    koverw       = Array{eltype(distr)}(undef, n_par.nk .* n_par.nm)
    for k = 1:n_par.nk
        for m = 1:n_par.nm
            total_wealth[m+(k-1)*n_par.nm] = n_par.grid_m[m] .+ Q .* n_par.grid_k[k]
            koverw[m+(k-1)*n_par.nm] = Q .* n_par.grid_k[k] ./ total_wealth[m+(k-1)*n_par.nm]
        end
    end
    # Wealth shares and gini
    IX = sortperm(total_wealth)
    total_wealth = total_wealth[IX]
    total_wealth_pdf = sum(distr, dims = 3)
    total_wealth_pdf = total_wealth_pdf[IX]
    total_wealth_cdf = cumsum(total_wealth_pdf)
    total_wealth_w = total_wealth .* total_wealth_pdf # weighted
    wealthshares = cumsum(total_wealth_w) ./ sum(total_wealth_w)
    shares = mylinearinterpolate(total_wealth_cdf, wealthshares, [0.5,0.9])
    BOTTOM50Wshare = shares[1]
    TOP10Wshare = 1.0 - shares[2]
    GiniW = gini(total_wealth, total_wealth_pdf)
    quantiles  = mylinearinterpolate(total_wealth_cdf,total_wealth.-minimum(total_wealth),[0.1,0.5,0.9])
    rat91W = quantiles[3]/quantiles[1]
    rat51W = quantiles[2]/quantiles[1]
    rat95W = quantiles[3]/quantiles[2]

    IX = sortperm(koverw)
    koverw = koverw[IX]
    koverw_pdf = sum(distr, dims = 3)
    koverw_pdf = koverw_pdf[IX]
    koverw_cdf = cumsum(koverw_pdf)
    koverw_w = koverw .* koverw_pdf # weighted
    KWshares = cumsum(koverw_w) ./ sum(koverw_w);
    shares = mylinearinterpolate(koverw_cdf,KWshares,[0.5,0.9])
    BOTTOM50KWshare = shares[1]
    TOP10KWshare = 1.0 - shares[2]

    capital_cdf = cumsum(distr_k)
    capital_w = n_par.grid_k .* distr_k # weighted
    Kshares = cumsum(capital_w) ./ sum(capital_w)
    shares = mylinearinterpolate(capital_cdf,Kshares,[0.5,0.9])
    BOTTOM50Kshare = shares[1]
    TOP10Kshare = 1.0 - shares[2]

    # Consumption distribution
    c = Array{eltype(c_a_star)}(undef, (n_par.nm, n_par.nk, n_par.ny, 2))
    distr_c = similar(c)
    aux_x = inc[5] # adjustment for labor in GHH preferences
    aux_x[:, :, end] .= zeros(n_par.nm, n_par.nk) # entrepreneurs do not work
    c[:, :, :, 1] .= c_a_star .+ aux_x # add adjustment to consumption
    c[:, :, :, 2] .= c_n_star .+ aux_x # add adjustment to consumption
    distr_c[:, :, :, 1] .= m_par.λ .* distr
    distr_c[:, :, :, 2] .= (1 - m_par.λ) .* distr

    # Gini of goods consumption
    IX = sortperm(c[:])
    c = c[IX]
    c_pdf = distr_c[IX]
    c_cdf = cumsum(c_pdf)
    Cshares = cumsum(c_pdf.*c)./sum(c.*c_pdf);
    shares = mylinearinterpolate(c_cdf, Cshares, [0.5,0.9])
    BOTTOM50Cshare = shares[1]
    TOP10Cshare = 1.0 - shares[2]
    GiniC = gini(c, c_pdf)
    quantiles = mylinearinterpolate(c_cdf,c[:],[0.1,0.5,0.9])
    rat91C = quantiles[3]/quantiles[1]
    rat51C = quantiles[2]/quantiles[1]
    rat95C = quantiles[3]/quantiles[2]

    # Top 10 net income share
    capital_inc = inc[2] .+ inc[3] .- n_par.mesh_m
    Yidio = inc[6] .+ capital_inc
    IX = sortperm(Yidio[:])
    Yidio = Yidio[IX]
    Y_pdf = distr[IX]
    Y_cdf = cumsum(Y_pdf)
    Y_w = Yidio .* Y_pdf
    net_incomeshares = cumsum(Y_w) ./ sum(Y_w)
    shares = mylinearinterpolate(Y_cdf, net_incomeshares, [0.5,0.9])
    BOTTOM50Inetshare = shares[1]
    TOP10Inetshare = 1.0 - shares[2]
    GiniInet = gini(Yidio, Y_pdf)
    quantiles  = mylinearinterpolate(Y_cdf,Yidio,[0.1,0.5,0.9])
    rat91Inet = quantiles[3]/quantiles[1]
    rat51Inet = quantiles[2]/quantiles[1]
    rat95Inet = quantiles[3]/quantiles[2]

    # Top 10 gross income share
    Yidio = incgross[1] .+ capital_inc
    IX = sortperm(Yidio[:])
    Yidio = Yidio[IX]
    Y_pdf = distr[IX]
    Y_cdf = cumsum(Y_pdf)
    Y_w = Yidio .* Y_pdf
    incomeshares = cumsum(Y_w) ./ sum(Y_w)
    shares = mylinearinterpolate(Y_cdf, incomeshares, [0.5,0.9])
    BOTTOM50Ishare = shares[1]
    TOP10Ishare = 1.0 - shares[2]
    GiniI = gini(Yidio, Y_pdf)
    quantiles  = mylinearinterpolate(Y_cdf,Yidio,[0.1,0.5,0.9])
    rat91I = quantiles[3]/quantiles[1]
    rat51I = quantiles[2]/quantiles[1]
    rat95I = quantiles[3]/quantiles[2]

    # Standard deviation of log labor earnings
    Yidio = log.(incgross[1][:, :, 1:end-1])
    IX = sortperm(Yidio[:])
    Yidio = Yidio[IX]
    distr_aux = distr[:, :, 1:end-1]
    distr_aux = distr_aux ./ sum(distr_aux[:])
    Y_pdf = distr_aux[IX]

    sdlogy = sqrt(dot(Y_pdf, Yidio .^ 2) .- dot(Y_pdf, Yidio) .^ 2)



    return distr_m, distr_k, distr_y,
        GiniW, TOP10Wshare, BOTTOM50Wshare, rat91W, rat51W, rat95W,
        GiniC, TOP10Cshare, BOTTOM50Cshare, rat91C, rat51C, rat95C,
        GiniInet, TOP10Inetshare, BOTTOM50Inetshare, rat91Inet, rat51Inet, rat95Inet, 
        GiniI, TOP10Ishare, BOTTOM50Ishare, rat91I, rat51I, rat95I,
        TOP10KWshare, BOTTOM50KWshare,
        TOP10Kshare, BOTTOM50Kshare,
        sdlogy
end
function gini(x, pdf)
    s = 0.0
    gini = 0.0
    for i in eachindex(x)
        gini -= pdf[i] * s
        s += x[i] * pdf[i]
        gini -= pdf[i] * s
    end
    gini /= s
    gini += 1.0
    return gini
end
