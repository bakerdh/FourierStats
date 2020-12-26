function importmousedata

% function to import mouse SSVEP data from the Hwang et al (2019) study
% data are stored in EEGlab format and have already been epoched
% the Fourier transform is taken, and complex spectra at the signal frequencies are stored for analysis in R
% DHB 21/12/20

close all;

datadir = '/Users/danbaker/Desktop/10.12751_g-node.e5tyek/';
nconds = [5 7];
startpoint = 1602;
endpoint = 3601;    % this is 1 second of data during which the stimulus was presented

mousedata = zeros(6,38,12,100);
for subj = 1:6
    for expt = 1:2
        
        EEG = pop_loadset(strcat('epochs_animal',num2str(subj),'.set'),strcat(datadir,'dataset_',num2str(expt),'/'));
        duration = EEG.pnts/EEG.srate;
        for ch = 1:EEG.nbchan
            elspec = [];
            for cond = 1:nconds(expt)
                
                condcount = 0;
                for n = 1:EEG.trials
                    thiscond = str2num(EEG.event(n).type(1));
                    if thiscond==cond
                        condcount = condcount + 1;
                        temp = fft(squeeze(EEG.data(ch,startpoint:endpoint,n)))/2000;
                        elspec(condcount,:) = temp(2:101);
                    end
                    
                end
                mousedata(subj,ch,(cond+(expt-1)*5),:) = mean(elspec,'omitnan');
            end
            
        end
    end
end


save('Hwangdata.mat','mousedata');

end