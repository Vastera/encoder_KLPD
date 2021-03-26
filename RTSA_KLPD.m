function [M_opt,kurt]=RTSA_KLPD(Phase,fs,fault_type,f_s,varargin)
% Copyright@ vastera@163.com
% General introduction: Find the optimal order of local polynomial diffferentiator, and record every corresponding kurtosis
%% ====================== INPUT ========================
% Phase:          Type:vector
%                          Phase description: input instantaneous phases
% fs:          Type: vector
%                          fs description: sampling frequency
% fault_type:         Type: string
%                    Freq description: fault type for choosingtime series for time synchronous averaging
% f_s:                Type: double
%                    f_s description: input sun gear rotating frequency
% ---------------------OPTIONAL:
%   filer_range:      Type: [LowerBound UpperBound] 1*2 vector
%                         description: the actual frequency range for filteration before envelope analysis
%   M_max:            Type:an odd number (default=401)
%                            description: the max length of iteration
%% ====================== OUTPUT =======================
% M_opt:          Type: an integer
%                           output_arg description: the optimal length of filter:L=2*M+1
% kurt:           Type: vector
%                       kurt description: all kurtosis result from the iteration of LPD for different L
%% =====================================================
narginchk(4, 8);
M_max = 401;
filter_range=[];
sig_env=true;
refs=4096;
[f_c,f_p_s,f_i,f_o,f_cg,f_b,f_m,f_sf,f_pf,f_rf] = CharacteristicFreq(f_s,36,35,108,3,0.0035,0.0195,10);
% Freq_blocked=[f_c:f_c:refs/2,2:2:refs/2,f_m:f_m:refs/2,f_s:f_s:refs/2];
if nargin>=5
    for argindex=1:2:length(varargin)
        if strcmp(varargin{argindex},'M_max')
            M_max =varargin{argindex+1};
        elseif strcmp(varargin{argindex},'filter_range')
            filter_range = varargin{argindex+1};
        elseif strcmpi(varargin{argindex},'sig_env')
            sig_env= varargin{argindex+1};
        elseif strcmpi(varargin{argindex},'Freq_blocked')
            Freq_blocked= varargin{argindex+1};
        end 
    end
end
T=length(Phase)/fs;
kurt=zeros(M_max,1);
h_wait=waitbar(0,'KLPD please wait...');
for m=1:M_max
    waitbar(m/M_max,h_wait,['KLPD please wait...',num2str(m),'/',num2str(M_max)])
    v=resample(LPD(Phase,m,fs),refs,fs);
    v=LinearPredict(v,100);
    if ~isempty(filter_range)
        b=fir1(256,filter_range/fs*2);
        v=filtfilt(b,1,v);
    end
    if sig_env
        v_env=envelope(v);
%         v_resam=resample(v_env,refs,fs);
        % f_s=abs(mean(v));
        %%%%%%%%%% BLock operation-related frequencies:rotation(f_s), meshing(f_m), carrier(f_c), frequency convertor(2) %%%%%%%
%         v_recon=real(Rm_Freq(v_resam,Freq_blocked,refs,2/T));% 3/T menas to block three points on the frequency axis each side.
    else
        v_resam=resample(v,refs,fs);
        v_recon=real(Rm_Freq(v_resam,Freq_blocked,refs,2/T));
    end
   
    %%%%%%%%%%%%%%%%% TSA %%%%%%%%%%%%%%%%%%%%%%%%%
     if strcmpi(fault_type,'sun')
        tp=0:2/f_sf:T;
    elseif strcmpi(fault_type,'ring')
        tp=0:2/f_rf:T;
    elseif strcmpi(fault_type,'planet')
        tp=0:2/f_pf:T;
    elseif strcmpi(fault_type,'inner')
        tp=0:2/f_i:T;
    elseif strcmpi(fault_type,'outer')
        tp=0:2/f_o:T;
    elseif strcmpi(fault_type,'RE')
        tp=0:2/f_b:T;
    end
    ta=tsa(v_env,refs,tp);
    kurt(m)=kurtosis(ta);
end
close(h_wait);
[~,M_opt]=max(kurt);
end