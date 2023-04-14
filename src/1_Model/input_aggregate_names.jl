# This file defines the sets of aggregate shocks, states (inluding shocks), and controls
# The code checks that the number of equations in the aggregate model is equal to the number 
# of aggregate variables excluding the distributional summary statistics. The latter are not 
# contained in the aggregate model code as they are parameter free but change whenever the 
# distribution changes and do not show up in any aggregate model equation.

shock_names = [:Z, :ZI, :μ, :μw, :A, :Rshock, :Gshock, :Tprogshock, :Sshock]


state_names = [
    "A",
    "Z",
    "ZI",
    "RB",
    "μ",
    "μw",
    "σ",
    "Ylag",
    "Bgovlag",
    "Tlag",
    "Ilag",
    "wlag",
    "qlag",
    "Clag",
    "av_tax_ratelag",
    "τproglag",
    "qΠlag",
    "Gshock",
    "Tprogshock",
    "Rshock",
    "Sshock",
]

# List cross-sectional controls / distributional summary variables (no equations in aggregate model expected)
distr_names = ["GiniW", "TOP10Wshare", "BOTTOM50Wshare", "rat91W", "rat51W", "rat95W", 
               "GiniC", "TOP10Cshare", "BOTTOM50Cshare", "rat91C", "rat51C", "rat95C",
               "GiniInet", "TOP10Inetshare", "BOTTOM50Inetshare", "rat91Inet", "rat51Inet", "rat95Inet", 
               "GiniI", "TOP10Ishare", "BOTTOM50Ishare", "rat91I", "rat51I", "rat95I",
               "TOP10KWshare", "BOTTOM50KWshare",
               "TOP10Kshare", "BOTTOM50Kshare",
               "sdlogy"]

control_names = [
    "r",
    "w",
    "K",
    "π",
    "πw",
    "Y",
    "C",
    "q",
    "N",
    "mc",
    "mcw",
    "u",
    "qΠ",
    "firm_profits",
    "RL",
    "Bgov",
    "Ht",
    "av_tax_rate",
    "T",
    "I",
    "B",
    "BD",
    "BY",
    "TY",
    "mcww",
    "G",
    "τlev",
    "τprog",
    "Ygrowth",
    "Bgovgrowth",
    "Igrowth",
    "wgrowth",
    "Cgrowth",
    "Tgrowth",
    "LP",
    "LPXA",
    "unionprofits",
    "profits",
]

# All controls in one array
control_names = [distr_names; control_names]
# All names in one array
aggr_names = [state_names; control_names]

# ascii names used for cases where unicode doesn't work, e.g., file saves
unicode2ascii(x) =
    replace.(
        replace.(
            replace.(replace.(replace.(x, "τ" => "tau"), "σ" => "sigma"), "π" => "pi"),
            "μ" => "mu",
        ),
        "ρ" => "rho",
    )

state_names_ascii = unicode2ascii(state_names)
control_names_ascii = unicode2ascii(control_names)
aggr_names_ascii = [state_names_ascii; control_names_ascii]
