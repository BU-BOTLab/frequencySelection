function ave_dat=averageASCData(fnames,numDiodes,stderr)
% fnames must have entire name and path
% Outputs: 
% ave_dat.error
% ave_dat.freq
% ave_dat.phase 
% ave_dat.AC
% ave_dat.dist
% ave_dat.ID
% ave_dat.diode_names
% ave_dat.nDiodes
% ave_dat.ACsd
% ave_dat.phsd
% Example
% ave_dat=averageASCData(fnames,diodes_selected,stderr)
% averageFDPMDataAtDiodes(char({'fdpmdata-1-dcswitch.asc'}),[658 830], [.03, .3])

diod(numDiodes) = struct('AC',[],'phase',[]);  %%array preallocation not really correct
for i=1:length(fnames)
    dat=getASCData(fnames{i});
    if dat.error~=0
        ave_dat.error=dat.error;
        return;
    end
    
    %Copy the interesting parts of the data files into the diod array
    for j = 1:dat.nDiodes
        diod(j).AC(:,i) = dat.AC(:,j);
        diod(j).phase(:,i) = dat.phase(:,j);
        diod(j).real(:,i) = dat.real(:,j);
        diod(j).imag(:,i) = dat.imag(:,j);
        try
            diod(j).posX = dat.posX;
            diod(j).posY = dat.posY;
        catch
            diode(j).posX = 0;
            diode(j).posY = 0;
        end
    end
    ave_dat.timestamp = dat.timestamp; %byh trying to avoid error with multiple files
    
    %Error checking
    if i==1 %If it's the first file just copy it
        ave_dat = dat;
    else
        if ave_dat.nDiodes ~= dat.nDiodes
            disp ('Program aborted:  inconsistent numbers of diodes between files to be averaged');
            ave_dat.error=-1;
            return;
        elseif  ave_dat.freq(1) ~= dat.freq(1)
            disp ('Program aborted:  inconsistent initial frequencies between files to be averaged');
            ave_dat.error=-1;
            return;
        elseif   length(ave_dat.freq) ~= length(dat.freq)
            disp ('Program aborted:  inconsistent numbers of frequencies between files to be averaged');
            ave_dat.error=-1;
            return;
        end     
    end   
end


for j = 1:ave_dat.nDiodes	
	%For each diode average the amplitude/phase/real/imaginary data
    ave_dat.AC(:,j) = mean(diod(j).AC,2);
	ave_dat.phase(:,j) = mean(diod(j).phase,2);
    ave_dat.real(:,j) = mean(diod(j).real,2);
    ave_dat.imag(:,j) = mean(diod(j).imag,2);
    %If there is more than one file, built up the error model of the
    %instrument
    if length(fnames)>1
        ave_dat.ACsd(:,j) = std(diod(j).AC,0,2);
		ave_dat.phsd(:,j) = std(diod(j).phase,0,2);
        ave_dat.realsd(:,j) = std(diod(j).real,0,2);
        ave_dat.imagsd(:,j) = std(diod(j).imag,0,2);
    end
end
%If there's only one file use the standard error values input into the
%function (usually 3% in amplitude and 0.3 degrees in phase)
if length(fnames)==1
    ave_dat.ACsd = stderr(1) .* ave_dat.AC;
	ave_dat.phsd = stderr(2) .*ones(size(ave_dat.phase));
    ave_dat.realsd = stderr(1) .* ave_dat.real;
    ave_dat.imagsd = stderr(1) .* ave_dat.imag;
    ave_dat.posX = dat.posX;
    ave_dat.posY = dat.posY;
else
    ave_dat.posX = mean([diod(:).posX]);
    ave_dat.posY = mean([diod(:).posY]);
end

%Convert from degrees to radians
ave_dat.phase=ave_dat.phase.*pi/180;  
ave_dat.phsd=ave_dat.phsd.*pi/180;   
ave_dat.error=0;