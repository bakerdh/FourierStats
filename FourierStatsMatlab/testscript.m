
load humanSSVEPdata.mat;

subjmeans = squeeze(alltarget(:,31,:));

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


data = subjmeans(:);
participantlabels = repmat(1:100,[7 1])';
participantlabels = participantlabels(:);
grouplabels = repmat(1:7,[100 1]);
grouplabels = grouplabels(:);

% anovacirc_test(data,grouplabels,participantlabels)
% anovacirc_test(data,grouplabels,[])

mergeddata = [real(data),imag(data),grouplabels,participantlabels];
% anovacirc_test(mergeddata,[],[])
mergeddata = [real(data),imag(data),grouplabels];
% anovacirc_test(mergeddata,[],[])

% analysecplx(subjmeans(:,1),[],[])

analysecplx(mergeddata,[],[])

% output = pairwisemahal(data,grouplabels)
% output.groups
% output.D

