function [hz_intermodulatory]=intermod(Hz_stimulation,limit)


    

% hz_stim=[3.5294    6.3158    8.0000];
% Hz_stimulation= [Hz_stimulation Hz_stimulation*2];

% index = randperm(length(Hz_stimulation));
% Hz_stimulation = Hz_stimulation(index);


ll = length(Hz_stimulation);
    hz_INTER = zeros(1,ll^3);                           % init hz_resp
    hz_INTER(1:ll) = (Hz_stimulation);                     % stimulation frequencies; imported to sort them
    if ll>1,
        kk = 0;
        for cc = 1:ll-2,
            for aa = 1:ll-1,
                for bb = aa+1:ll,
                    kk = kk + 1;
                    hz_INTER(ll+2*(kk-1)+1) = abs(hz_INTER(bb)-hz_INTER(aa)); % intermodulary frequencies; uses the sorted stimulation frequencies
                    hz_INTER(ll+2*(kk-1)+2) = abs(hz_INTER(bb)+hz_INTER(aa));
                    hz_INTER(ll+3*(kk-1)+3) = abs(hz_INTER(bb)+hz_INTER(aa)+hz_INTER(cc));
%                     hz_INTER(ll+4*(kk-1)+4) = abs(hz_INTER(bb)+hz_INTER(aa)-hz_INTER(cc));
                    hz_INTER(ll+5*(kk-1)+5) = abs(2*hz_INTER(bb)-hz_INTER(aa));
                    hz_INTER(ll+6*(kk-1)+6) = abs(2*hz_INTER(bb)-2*hz_INTER(aa));
                    
%                     hz_INTER(ll+5*(kk-1)+5) = abs(hz_INTER(bb)+2*hz_INTER(aa));
%                     hz_INTER(ll+3*(kk-1)+3) = abs(hz_INTER(bb)-hz_INTER(aa)+hz_INTER(cc));
%                     hz_INTER(ll+3*(kk-1)+3) = abs(hz_INTER(bb)+hz_INTER(aa)+hz_INTER(cc));
                end
            end
        end
    end
    
  hz_INTER(hz_INTER==0)=[];
  hz_INTER(hz_INTER>limit(2))=[];
  hz_INTER(hz_INTER<limit(1))=[];
  hz_INTER=unique(hz_INTER,'stable');
  idx = ismember(hz_INTER,Hz_stimulation);
  hz_intermodulatory=hz_INTER(~idx);
    
    
end