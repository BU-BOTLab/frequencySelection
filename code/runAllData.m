%runs all of the scripts to generate all of the data in the paper
addpath(genpath('utilities'));
fprintf("Generating OP test set\n");
getOPSet;
%Noise free simulations
fprintf("Working on noise free simulations\n")
noNoise_LMA_default;
noNoise_LMA_optimized;
noNoise_LUT;

%All-digital noise model
fprintf("Working on the all-digital system simulations\n")
fprintf("\tEqual data\n")
fprintf("\tEqual data, default\n")
Digital_noisemodel_equalData_LMA_default;
fprintf("\tEqual data, optimized\n")
Digital_noisemodel_equalData_LMA_optimized;
fprintf("\tEqual data, LUT\n")
Digital_noisemodel_equalData_LUT;
fprintf("\tUnequal data, default\n")
Digital_noisemodel_unequalData_LMA_default;
fprintf("\tUnequal data, optimized\n")
Digital_noisemodel_unequalData_LMA_optimized;
fprintf("\tUnequal data, LUT\n")
Digital_noisemodel_unequalData_LUT;
Digital_noisemodel_equalData_HDLUT_singleFreq;

%Network analyzer noise model
fprintf("Working on the network analyzer system simulations\n")
fprintf("\tEqual data\n")
NA_noisemodel_equalData_LMA_default;
NA_noisemodel_equalData_LMA_optimized;
NA_noisemodel_equalData_LUT;
fprintf("\tUnequal data\n")
NA_noisemodel_unequalData_LMA_default;
NA_noisemodel_unequalData_LMA_optimized;
NA_noisemodel_unequalData_LUT;

%Other experimental runs
NA_noisemodel_equalData_HDLUT_singleFreqs;
compareInversionTiming;
LUTerrorVsDensity_1freq;
findOptimalFreq;
getSDVsOP;