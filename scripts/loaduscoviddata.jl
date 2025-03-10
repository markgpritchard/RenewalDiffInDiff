# Data used by Wei Lyu and George L. Wehby, https://doi.org/10.1377/hlthaff.2020.00818

casesdf = CSV.read(
    datadir("exp_raw", "us-states.csv"), DataFrame
)
# remove data after 22 May 
filter!(:date => x -> x <= Date("2020-05-22"), casesdf)
# make a DataFrame of states and fips identities

statepopulationsdf = CSV.read(
    datadir(
        "exp_raw", 
        "COVID-19 US state policy database (CUSP) - State Characteristics.csv"
    ), DataFrame;
    footerskip=1  # final row sums some variables
)

facemasksdf = CSV.read(
    datadir(
        "exp_raw", 
        "COVID-19 US state policy database (CUSP) - Face Masks.csv"
    ), DataFrame;
    footerskip=1  # final row sums some variables
)

# different numbers of states 
filter!(:fips => x -> x ∈ getproperty(
    statepopulationsdf, Symbol("State FIPS Code")), casesdf
)

states = select(casesdf, :state, :fips)
unique!(states)
sort!(states)

incidence = zeros(Int, 123, 51)
for (i, state) ∈ enumerate(states.state)
    stateinds = findall(x -> x == state, casesdf.state)
    for (j, d) ∈ enumerate(Date("2020-01-21"):Date("2020-05-22"))
        rawdateind = findfirst(x -> x == d, casesdf.date[stateinds])
        isnothing(rawdateind) && continue
        incidence[j, i] = (
            casesdf.cases[stateinds[rawdateind]] - 
            (j == 1 ? 0 : incidence[j-1, i])
        ) 
    end
end

populations = getproperty(statepopulationsdf, Symbol("Population 2018"))

W_uscoviddata = generatew_gt(COVIDSERIALINTERVAL, incidence, populations)

maskday = (
    datestrings = getproperty(facemasksdf, Symbol("Public face mask mandate start"));
    dates = [
        x == "0" ? 
            Date("2099-12-31") : 
            Date(x, "m/d/yyyy")
        for x ∈ datestrings 
    ];
    starttimeperiods = dates .- Date("2020-01-21");
    starttimes = Dates.value.(starttimeperiods);
    InterventionsMatrix{Int}(starttimes, 123)
)

relaxshelterinplace = (
    df = CSV.read(
        datadir("exp_raw", "COVID-19 US state policy database (CUSP) - Stay at Home.csv"), 
        DataFrame;
        footerskip=1,
    );
    datestrings = getproperty(df, Symbol("End stay at home/shelter in place "));
    dates = [
        x == "0" ? 
            Date("2099-12-31") : 
            Date(x, "m/d/yyyy")
        for x ∈ datestrings 
    ];
    starttimeperiods = dates .- Date("2020-01-21");
    starttimes = Dates.value.(starttimeperiods);
    InterventionsMatrix{Int}(starttimes, 123)
)

_openandclosedf = CSV.read(
    datadir(
        "exp_raw", 
        "COVID-19 US state policy database (CUSP) - Closures & Reopening.csv"
    ), 
    DataFrame;
    footerskip=1,
);

reopenbusiness = (
    datestrings = getproperty(
        _openandclosedf, 
        Symbol("Began to reopen businesses statewide")
    );
    dates = [
        x == "0" ? 
            Date("2099-12-31") : 
            Date(x, "m/d/yyyy")
        for x ∈ datestrings 
    ];
    starttimeperiods = dates .- Date("2020-01-21");
    starttimes = Dates.value.(starttimeperiods);
    InterventionsMatrix{Int}(starttimes, 123)
)

reopenrestaurants = (
    datestrings = getproperty(
        _openandclosedf, 
        Symbol("Reopened restaurants")
    );
    dates = [
        x == "0" ? 
            Date("2099-12-31") : 
            Date(x, "m/d/yyyy")
        for x ∈ datestrings 
    ];
    starttimeperiods = dates .- Date("2020-01-21");
    starttimes = Dates.value.(starttimeperiods);
    InterventionsMatrix{Int}(starttimes, 123)
)

reopengyms = (
    datestrings = getproperty(
        _openandclosedf, 
        Symbol("Reopened gyms")
    );
    dates = [
        x == "0" ? 
            Date("2099-12-31") : 
            Date(x, "m/d/yyyy")
        for x ∈ datestrings 
    ];
    starttimeperiods = dates .- Date("2020-01-21");
    starttimes = Dates.value.(starttimeperiods);
    InterventionsMatrix{Int}(starttimes, 123)
)

reopencinemas = (
    datestrings = getproperty(
        _openandclosedf, 
        Symbol("Reopened movie theaters")
    );
    dates = [
        x == "0" ? 
            Date("2099-12-31") : 
            Date(x, "m/d/yyyy")
        for x ∈ datestrings 
    ];
    starttimeperiods = dates .- Date("2020-01-21");
    starttimes = Dates.value.(starttimeperiods);
    InterventionsMatrix{Int}(starttimes, 123)
)
