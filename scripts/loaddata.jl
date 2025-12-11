

using DrWatson
@quickactivate :RenewalDiffInDiff

import CSV 
using DataFrames
using Dates
using RenewalDiD

# UK masking data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

maskcoviddf = loadukmaskdata(
    datadir("exp_raw", "OxCGRT_GBR_differentiated_withnotes_2020.csv")
)
observedcases = ukobservedcasesmatrix(maskcoviddf)

# do not keep fitting parameters after all groups are treated: 
#   day 160 for masksrecommended
#   day 227 for masksrequired

# H6E_FacialCoverings == 1 means recommended 
masksrecommended = ukinterventions(maskcoviddf, 160, :H6E_FacialCoverings, 1)  
# H6E_FacialCoverings >= 2 means some level of requirement
masksrequired = ukinterventions(maskcoviddf, 226, :H6E_FacialCoverings, 2) 
# C1E_Schoolclosing >= 2 means require some school closing, 1 means recommend, 0 means no 
# measures
# no groups have school reopening before day 160
schoolreopening_226 = ukinterventions(maskcoviddf, 226, :C1E_Schoolclosing, x -> x <= 1, 2) 
# C6E_Stayathome == 0 means no stay at home measures
endstayathome_160 = ukinterventions(maskcoviddf, 160, :C6E_Stayathome, x -> x == 0, 2)    
endstayathome_226 = ukinterventions(maskcoviddf, 226, :C6E_Stayathome, x -> x == 0, 2)    
# H3E_Contacttracing ==2 means "comprehensive contact tracing; done for all identified cases"
comprehensivecontacttracing_160 = ukinterventions(maskcoviddf, 160, :H3E_Contacttracing, 2)    
comprehensivecontacttracing_226 = ukinterventions(maskcoviddf, 226, :H3E_Contacttracing, 2)    

# do not keep fitting parameters after all groups are treated: 
#   day 160 for masksrecommended
#   day 227 for masksrequired

covidmaskdata1 = RenewalDiDData( ;
    observedcases=observedcases[1:161, :],
    interventions=masksrecommended,
    Ns=UKPOPULATION2020,
    minvalue=0.1, 
    sampletime=10,
    id="data1: masks recommended, no other interventions"
)
safesave(datadir("exp_pro", "covidmaskdata1.jld2"), Dict("data" => covidmaskdata1))

covidmaskdata2 = RenewalDiDData( ;
    observedcases=observedcases[1:227, :],
    interventions=masksrequired,
    Ns=UKPOPULATION2020,
    minvalue=0.1, 
    sampletime=10,
    id="data2: masks required, no other interventions"
)
safesave(datadir("exp_pro", "covidmaskdata2.jld2"), Dict("data" => covidmaskdata2))

covidmaskdata3 = RenewalDiDData( ;
    observedcases=observedcases[1:161, :],
    interventions=InterventionArray(masksrecommended; offset=-35:7:35),
    Ns=UKPOPULATION2020,
    minvalue=0.1, 
    sampletime=10,
    id="data3: masks recommended, lead and lag placebos"
)
safesave(datadir("exp_pro", "covidmaskdata3.jld2"), Dict("data" => covidmaskdata3))

covidmaskdata4 = RenewalDiDData( ;
    observedcases=observedcases[1:227, :],
    interventions=InterventionArray(masksrequired; offset=-35:7:35),
    Ns=UKPOPULATION2020,
    minvalue=0.1, 
    sampletime=10,
    id="data4: masks required, lead and lag placebos"
)
safesave(datadir("exp_pro", "covidmaskdata4.jld2"), Dict("data" => covidmaskdata4))

covidmaskdata5 = RenewalDiDData( ;
    observedcases=observedcases[1:161, :],
    interventions=cat(
        masksrecommended, endstayathome_160, comprehensivecontacttracing_160; 
        dims=3
    ),
    Ns=UKPOPULATION2020,
    minvalue=0.1, 
    sampletime=10,
    id="data5: masks recommended, competing interventions"
)
safesave(datadir("exp_pro", "covidmaskdata5.jld2"), Dict("data" => covidmaskdata5))

covidmaskdata6 = RenewalDiDData( ;
    observedcases=observedcases[1:227, :],
    interventions=cat(
        masksrequired,
        schoolreopening_226, 
        endstayathome_226, 
        comprehensivecontacttracing_226; 
        dims=3
    ),
    Ns=UKPOPULATION2020,
    minvalue=0.1, 
    sampletime=10,
    id="data6: masks required, competing interventions"
)
safesave(datadir("exp_pro", "covidmaskdata6.jld2"), Dict("data" => covidmaskdata6))

covidmaskdata7 = RenewalDiDData( ;
    observedcases=observedcases[1:161, :],
    interventions=(
        iv = InterventionArray(masksrecommended; offset=-35:7:35);
        cat(iv, endstayathome_160, comprehensivecontacttracing_160; dims=3)
    ),
    Ns=UKPOPULATION2020,
    minvalue=0.1, 
    sampletime=10,
    id="data7: masks recommended, competing interventions and lead and lag placebos"
)
safesave(datadir("exp_pro", "covidmaskdata7.jld2"), Dict("data" => covidmaskdata7))

covidmaskdata8 = RenewalDiDData( ;
    observedcases=observedcases[1:227, :],
    interventions=(
        iv = InterventionArray(masksrequired; offset=-35:7:35);
        cat(
            iv, schoolreopening_226, endstayathome_226, comprehensivecontacttracing_226; 
            dims=3
        )
    ),
    Ns=UKPOPULATION2020,
    minvalue=0.1, 
    sampletime=10,
    id="data8: masks required, competing interventions and lead and lag placebos"
)
safesave(datadir("exp_pro", "covidmaskdata8.jld2"), Dict("data" => covidmaskdata8))


