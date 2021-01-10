
load humanSSVEPdata.mat;

% subjmeans = squeeze(alltarget(:,31,:));

% for n = 1:7
%     temp = tsqc_test(subjmeans(:,n),[],[],[]);
%     alloutput(n) = temp.tsqc;
% end
% temp
% alloutput
% 
% 
% for n = 1:7
%     temp = tsqc_test(subjmeans(:,1),subjmeans(:,n),0,[]);
%     alloutput(n) = temp.tsqc;
% end
% temp
% alloutput
% 
% 
% for n = 1:7
%     temp = tsqc_test(subjmeans(:,1),subjmeans(:,n),1,[]);
%     alloutput(n) = temp.tsqc;
% end
% temp
% alloutput
% 
% 
% 
% for n = 1:7
%     temp = tsqh_test(subjmeans(:,n),[],[],[]);
%     alloutput(n) = temp.tsq;
% end
% temp
% alloutput
% 
% 
% for n = 1:7
%     temp = tsqh_test(subjmeans(:,1),subjmeans(:,n),0,[]);
%     alloutput(n) = temp.tsq;
% end
% temp
% alloutput
% 
% 
% for n = 2:7
%     temp = tsqh_test(subjmeans(:,1),subjmeans(:,n),1,[]);
%     alloutput(n) = temp.tsq;
% end
% temp
% alloutput
% 
% 
% 
% 
% for n = 1:7
%     xydata = [real(subjmeans(:,n)) imag(subjmeans(:,n))]';
%     temp = tsqc_test(xydata,[],[],[]);
%     alloutput(n) = temp.tsqc;
% end
% temp
% alloutput
% 
% 
% for n = 1:7
%     xydata = [real(subjmeans(:,n)) imag(subjmeans(:,n))]';
%     temp = tsqh_test(xydata,[],[],[]);
%     alloutput(n) = temp.tsq;
% end
% temp
% alloutput

% for n = 1:7
%     temp = CI_test(subjmeans(:,n),0.05);
%     alloutput(n) = temp.CI;
% end
% temp
% alloutput


% 
% for n = 1:7
%     xydata = [real(subjmeans(:,n)) imag(subjmeans(:,n))]';
%     temp = CI_test(xydata,0.05);
%     alloutput(n) = temp.CI;
% end
% temp
% alloutput


% data = subjmeans(:);
% participantlabels = repmat(1:100,[7 1])';
% participantlabels = participantlabels(:);
% grouplabels = repmat(1:7,[100 1]);
% grouplabels = grouplabels(:);

% anovacirc_test(data,grouplabels,participantlabels)
% anovacirc_test(data,grouplabels,[])

% mergeddata = [real(data),imag(data),grouplabels,participantlabels];
% anovacirc_test(mergeddata,[],[])
% mergeddata = [real(data),imag(data),grouplabels];
% anovacirc_test(mergeddata,[],[])

% analysecplx(subjmeans(:,1),[],[])

% analysecplx(mergeddata,[],[])

% output = pairwisemahal(data,grouplabels)
% output.groups
% output.D

        threshdist = 0.18;          % arbitrary distance between electrodes that seems to work well in pairing up neighbouring electrodes
        
        for m = 1:64                % loop through all electrodes
            for n = 1:64            % calculate absolute distance from target electrode to each other electrode
                xy1 = montage.electrodelocs(montage.channelmappings(m),:);
                xy2 = montage.electrodelocs(montage.channelmappings(n),:);
                electrodedistances(m,n) = sqrt((xy1(1)-xy2(1)).^2 + (xy1(2)-xy2(2)).^2);
            end
        end
        
        adjacencymatrix = zeros(size(electrodedistances));      
        adjacencymatrix(find(electrodedistances<threshdist)) = 1;   % electrode pairs that are less than the threshold distance away get set to 1
        adjacencymatrix(find(electrodedistances==0)) = 0;           % an electrode can't be paired with itself (set the major diagonal to 0)


% output = clustercorrect(squeeze(alltarget(:,1:64,1)),squeeze(alltarget(:,1:64,7)),adjacencymatrix,3);
% output = clustercorrect(squeeze(alltarget(:,1:64,4)),[],adjacencymatrix,3);
% a = output.clusterpoints{1}
% b = output.clusterpoints{2}
% output

datax = squeeze(alltarget(:,1:64,4));
clustformthresh = 0.00001;
output = clustercorrect(squeeze(alltarget(:,1:64,4)),[],adjacencymatrix,3,0,0.00001);
output
output.clusterpoints{1}
output.clusterpoints{2}

% data = randn(20,100);
% data(:,50:55) = data(:,50:55) + 0.5;
% data(:,20:25) = data(:,20:25) + 0.5;
% output = clustercorrect(data,[],[],1);
% 



