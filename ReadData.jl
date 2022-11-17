using JLD2
using Plots
using DelimitedFiles

function ReadData(fileName)
    if ContainsMeanAndVariance(fileName) && !(ContainsProperFluc(fileName))
        data = load("./Data for the paper/" * fileName * "/data.jld2", "x", "y")
        return data[1], SeparateMeanAndVar(data[2][1]), SeparateMeanAndVar(data[2][2]), SeparateMeanAndVar(data[2][3])
    else
        data = load("./Data for the paper/" * fileName * "/data.jld2", "x", "y")
        return data[1], data[2][1], data[2][2], data[2][3]
    end
end
function ReadData(foldername, fileName)
    if ContainsMeanAndVariance(fileName) && !(ContainsProperFluc(fileName))
        data = load("./" * foldername * "/" * fileName * "/data.jld2", "x", "y")
        return data[1], SeparateMeanAndVar(data[2][1]), SeparateMeanAndVar(data[2][2]), SeparateMeanAndVar(data[2][3])
    else
        data = load("./" * foldername * "/" * fileName * "/data.jld2", "x", "y")
        return data[1], data[2][1], data[2][2], data[2][3]
    end
end

function ContainsMeanAndVariance(fileName)
    fileNamesWithVector = ["halfLN", "LN", "smallLN", "corrlN"]
    return !(chop(fileName, head = 3, tail = 0) in fileNamesWithVector)
end

function ContainsProperFluc(filename)
    return chop(filename, head = 3, tail = 0) == "Proper_fluctuation"
end

function SeparateMeanAndVar(data)
    @show data[1]
    return [res[1] for res in data], [res[2] for res in data]
end

function SaveForMatlab(fileName)
    x, y4, y6, y8 = ReadData(fileName)
    SaveVectorMatlab(x, fileName, "probabilities")
    if ContainsProperFluc(fileName)
        SaveVectorMatlab(y12, fileName, "L12")
        SaveVectorMatlab(y4, fileName, "L4")
        SaveVectorMatlab(y6, fileName, "L6")
        SaveVectorMatlab(y8, fileName, "L8")
    elseif ContainsMeanAndVariance(fileName)
        SaveVectorMatlab(y4[1], fileName, "L4_Mean")
        SaveVectorMatlab(y4[2], fileName, "L4_Var")
        SaveVectorMatlab(y6[1], fileName, "L6_Mean")
        SaveVectorMatlab(y6[2], fileName, "L6_Var")
        SaveVectorMatlab(y8[1], fileName, "L8_Mean")
        SaveVectorMatlab(y8[2], fileName, "L8_Var")
    else
        SaveVectorMatlab(reduce(hcat, y4), fileName, "L4")
        SaveVectorMatlab(reduce(hcat, y6), fileName, "L6")
        SaveVectorMatlab(reduce(hcat, y8), fileName, "L8")
    end
end

function SaveVectorMatlab(data, observableName, dataTypeName)
    path = pwd() * "/Data Matlab/" * observableName
	mkpath(path)
    writedlm(path * "/" * dataTypeName * ".csv", data)
end

function f()
    #Example how to read the data from the .jld2 files
    #For stuff with mean and variance y[1] is the mean, y[2] is the variance
    #x, y12 = ReadData("M1_Entanglement")
    #display(length(y8[1]))
    #plot!(x, y12[1])
    #To save the data in .csv form for matlab
    #SaveForMatlab("M2_smallLN")
end

function combinedata()
    observables = ["corrlN", "corrlNexpval", "Entanglement", "fluctuations", "halfLN", "halfLNexpval", "LN", "LNexpval", "Proper_fluctuation", "smallLN", "smallLNexpval"]
    for msr in ["M2", "M3"]
        for obs in observables
            filename = msr * "_" * obs
            x, y1 = ReadData("data1", filename)
            x, y2 = ReadData("data2", filename)
            y = vcat(y1, y2)
            mkpath("./" * "newdata" * "/" * filename)
            jldsave("./" * "newdata" * "/" * filename * "/data.jld2"; x, y)
        end
    end
end
