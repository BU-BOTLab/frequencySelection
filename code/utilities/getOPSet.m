%getOPSet.m
%
%PURPOSE: Generate a set of 10k optical properties randomly distributed
%within the range of interest. Saves them as a mat file
%
%INPUTS: 
%
% muaAll:    Range of mua values
% musAll:    Range of mus values
% numPairs:  Number of OP pairs to generate
% precision: How many significant digits to round to 
%
%OUTPUTS:
%
% trueMua: random values of absorption
% trueMus: random values of scattering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng(14850); %Set random seed for repeatability
%Range of OPs
muaAll = [0.001,0.05];
musAll = [.1,3];
numPairs =10000;
precision = 6;
trueMua = zeros(1,numPairs);
trueMus = zeros(1,numPairs);
for i = 1:numPairs
    %Get random OPs
    randMua = round(muaAll(1) + (muaAll(2) - muaAll(1)) * rand(1),precision);
    randMus = round(musAll(1) + (musAll(2) - musAll(1)) * rand(1),precision);
    %Make sure they satisfy diffusion condition
    while (10*randMua) > randMus
        randMua = round(muaAll(1) + (muaAll(2) - muaAll(1)) * rand(1),precision);
        randMus = round(musAll(1) + (musAll(2) - musAll(1)) * rand(1),precision);
    end
    %Save the true OP pair
    trueMua(i) = randMua;
    trueMus(i) = randMus;
end

save('../generatedData/OPSet10k.mat','trueMua','trueMus')