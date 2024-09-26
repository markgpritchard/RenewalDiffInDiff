
betafunction1_1(t) = 0.6 * (t < 100 ? 1.0 : 1.1)
betafunction1_2_counterfactual(t) = 0.7 * (t < 100 ? 1.0 : 1.1)
betafunction1_2(t) = (t < 100 ? 1.0 : 0.8) * betafunction1_2_counterfactual(t)
betafunctions1 = [ betafunction1_1, betafunction1_2 ]
betafunctions1_counterfactual = [ betafunction1_1, betafunction1_2_counterfactual ]

betafunction2_1(t) = 0.9 * (0.6 + 0.2 * cos(2π * (t - 80) / 365))
betafunction2_2_counterfactual(t) = 0.8 * (0.6 + 0.2 * cos(2π * (t - 80) / 365))
betafunction2_2(t) = (t < 75 ? 1.0 : 0.8) * betafunction2_2_counterfactual(t)
betafunction2_3_counterfactual(t) = 1.1 * (0.6 + 0.2 * cos(2π * (t - 80) / 365))
betafunction2_3(t) = (t < 100 ? 1.0 : 0.8) * betafunction2_3_counterfactual(t)
betafunctions2 = [ betafunction2_1, betafunction2_2, betafunction2_3 ]
betafunctions2_counterfactual = [ 
    betafunction2_1, betafunction2_2_counterfactual, betafunction2_3_counterfactual 
]

betafunction3_1(t) = 1 * (0.6 + 0.2 * cos(2π * (t + 80) / 365))
betafunction3_2(t) = (t < 110 ? 1.0 : 1.5) * 0.9 * (0.6 + 0.2 * cos(2π * (t + 80) / 365))
betafunction3_3_counterfactual(t) = 1.2 * (0.6 + 0.2 * cos(2π * (t + 80) / 365))
betafunction3_3(t) = (t < 75 ? 1.0 : 0.8) * betafunction3_3_counterfactual(t)

function betafunction3_4_counterfactual(t)
    (t < 50 ? 1.0 : 1.5) * 0.8 * (0.6 + 0.2 * cos(2π * (t + 80) / 365))
end

betafunction3_4(t) = (t < 125 ? 1.0 : 0.8) * betafunction3_4_counterfactual(t)
betafunctions3 = [ betafunction3_1, betafunction3_2, betafunction3_3, betafunction3_4 ]
betafunctions3_counterfactual = [ 
    betafunction3_1, betafunction3_2, 
    betafunction3_3_counterfactual, betafunction3_4_counterfactual
]

betafunction4_1(t) = 0.9 * (0.6 + 0.2 * cos(2π * (t - 80) / 365))
betafunction4_2_counterfactual(t) = (0.8 * (0.6 + 0.2 * cos(2π * (t - 80) / 365)))^1.3
betafunction4_2(t) = (t < 75 ? 1.0 : 0.8) * betafunction4_2_counterfactual(t)
betafunction4_3_counterfactual(t) = (1.1 * (0.6 + 0.2 * cos(2π * (t - 80) / 365)))^1.5
betafunction4_3(t) = (t < 100 ? 1.0 : 0.8) * betafunction4_3_counterfactual(t)
betafunctions4 = [ betafunction4_1, betafunction4_2, betafunction4_3 ]
betafunctions4_counterfactual = [ 
    betafunction4_1, betafunction4_2_counterfactual, betafunction4_3_counterfactual 
]
