function [out,mua,mus] = read_dbMU(fname)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
fid=fopen(fname,'r');
out = struct('name','','wv',0','mua',0,'dmua',0,'mus',0,'dmus',0);
header = fgetl(fid);
line= fgetl(fid);
iter = 1;
while ischar(line)
    splitsies = strsplit(line,'\t');
    out(iter).name = splitsies{3};
    out(iter).wv = str2double(splitsies{4});
    out(iter).mua = str2double(splitsies{5});
    out(iter).dmua = str2double(splitsies{6});
    out(iter).mus = str2double(splitsies{7});
    out(iter).dmus = str2double(splitsies{8});
    iter=iter+1;
    line = fgetl(fid);
end
mua = zeros(1,length(out));
mus = zeros(1,length(out));
for i = 1:length(out)
    mua(i) = out(i).mua;
    mus(i) = out(i).mus;
end
fclose(fid);
