% Reads fdpm data files 
% This function should only be called by averageASCData, but you can use it
% if you want to. Just be careful of the units for phase

% Last modified: Sept. 5 2019 MBA
% Filename should include entire path
% Outputs: 
% dat.error
% dat.freq
% dat.phase 
% dat.real
% dat.imag
% dat.AC
% dat.dist
% dat.ID
% dat.diode_names
% dat.nDiodes
% data.timestamp
% Example Call: 
% getASCData('fdpmdata-1-dcswitch.asc')


function dat=getASCData(fname)

fid=fopen(fname,'r'); %Open the ASC file

if fid==-1
    fprintf('Could not open file %s.',fname);
    dat.error=-1;
    return
end

%Read the header and extract a bunch of info
fw='';
dat.posX = 0;
dat.posY = 0;
while ~strcmp(fw,'Frequency')   %this section reads needed header info
    l=fgetl(fid);
    fw=strtok(l);
    
    if strfind(l,'Source-Detector (mm):')
        dat.dist = str2double(strtok(l(22:length(l))));
    elseif strfind(l,'Src X')
        s1 = strsplit(l,', ');
        s2_1 = strsplit(s1{1}, ' ');
        s2_2 = strsplit(s1{2}, ' ');
        dat.srcX = str2double(s2_1{3});
        dat.srcY = str2double(s2_2{3});
    elseif strfind(l,'DetX')
        s1 = strsplit(l,', ');
        s2_1 = strsplit(s1{1}, ' ');
        s2_2 = strsplit(s1{2}, ' ');
        dat.detX = str2double(s2_1{2});
        dat.detY = str2double(s2_2{2});
    elseif strfind(l,'Patient ID:')
        dat.ID = strtok(l(12:length(l)));
    elseif strfind(l,'Laser names:')
        dat.wavelengths = sscanf(l(13:length(l)), '%d');
        if isempty(dat.wavelengths)
            splitsies=strsplit(l,',');
            for i = 2:length(splitsies)
                dat.wavelengths = [dat.wavelengths,str2double(splitsies{i}(1:3))];
            end
        end
        if length(dat.wavelengths) == 1 && dat.wavelengths > 1000
            dat.wavelengths = [658,690,785,808,830,850];
        end
    elseif strfind(l,'Number of Lasers')
        s = strsplit(l,':');
        dat.nDiodes = str2double(s{2});
    elseif strfind(l,'Date acquired:')
        dat.timestamp = strtok(l(26:length(l)));
    elseif strfind(l, 'Pos X: ')
        s = strsplit(l,':');
        k = strsplit(s{2},',');
        dat.posY = str2double(s{3});
        dat.posX = str2double(k{1});
    end
end
            
A = fscanf(fid,'%lg',[1+dat.nDiodes*2 inf]);   % read numerical data

fclose(fid);  

dat.freq=A(1,:)'; %Frequency
badRows = [];
for j=1:dat.nDiodes
	dat.phase(:,j) = A(2+ 2*(j-1), :)'; %Phase (degrees)
	dat.AC(:,j) = A(3+ 2*(j-1), :)'; %Amplitude
    br = find(isnan(dat.phase(:,j)));
    if ~isempty(br)
        badRows = [badRows, br];
    end
end
dat.phase(badRows,:) = [];
dat.AC(badRows,:) = [];
dat.freq(badRows,:) = [];
%dat.phase = dat.phase .* pi/180;
dat.real=dat.AC .* cos(dat.phase);
dat.imag = dat.AC .* sin(dat.phase);
dat.error=0;
end


