% This script reproduces the simulation figures from the paper:
% "A unified framework for dissecting the effects of common signals on
% functional and effective connectivity analyses: power, coherence, and
% Granger causality" by Dror Cohen and Naotsugu Tsuchiya, 2017.

% The script requires the Fieldtrip toolbox (http://www.fieldtriptoolbox.org/)
% to be included in the matlab path

% The code was tested on Matlab version 2016b and Fieldtrip version June
% 2017  (revision 4cb126f)

% please email dror.cohen07@gmail.com regarding any questions or comments

clear all;
close all;

% fix seed
rng(0)
clf
analytic_data=[];
% ar model
analytic_data.coeffs=[0.1 0.4;
                      0 0.1;]; %off-diagonal entries simulate 2->1
analytic_data.noisecov=[1 0;
                         0 1;];   
analytic_data.fsampleorig=400;
analytic_data.dimord='chan_chan_lag';
sigs={};
for i=1:size(analytic_data.coeffs,1)
    sigs{i}=num2str(i);
end
analytic_data.label=sigs;
 
analytic_cfg        = [];
analytic_cfg.method = 'mvar';
connected_mfreq      = ft_freqanalysis(analytic_cfg, analytic_data);

% take mean as the power of the common input
UU=squeeze(mean(connected_mfreq.crsspctrm(1,1,:),3)+mean(connected_mfreq.crsspctrm(2,2,:),3))/2;

% simulate data
cfg = [];
cfg.method = 'ar';
cfg.ntrials = 800;
cfg.triallength = 1;
cfg.fsample = 400;
cfg.nsignal = 2;
cfg.params=analytic_data.coeffs;
cfg.noisecov=analytic_data.noisecov; 
sim_data = ft_connectivitysimulation(cfg);

% frequency analysis of sim data
cfg1 = [];
cfg1.method = 'mtmfft';
cfg1.taper = 'dpss';
cfg1.output = 'fourier';
cfg1.tapsmofrq = 10;
empirical_freq = ft_freqanalysis(cfg1, sim_data);

% sim data with common input. This cane be used for comparisons against analytic expressions
sim_data_common_input=sim_data;
for trialIndx=1:length(sim_data.trial)
      this_trial_ref=sqrt(UU/2)*randn(1,size(sim_data_common_input.trial{trialIndx},2));
      this_trial_ref=repmat(this_trial_ref,size(sim_data_common_input.trial{trialIndx},1),1);
      sim_data_common_input.trial{trialIndx}=sim_data_common_input.trial{trialIndx}+this_trial_ref;
end

% freq analysis
empirical_freq_common_input = ft_freqanalysis(cfg1, sim_data_common_input);

%%%%%%%%%%%%%%%%%%%%%%%%
%% create the scenarios
% Scneario 1 - unidirectionally connected system
mfreq_scenarios = {};
mfreq_scenarios{1} = connected_mfreq;

% Scenario 2 - disconneceted system
tmp_mfreq = connected_mfreq;
tmp_mfreq.crsspctrm(1,2,:) = 0;
tmp_mfreq.crsspctrm(2,1,:) = 0;
mfreq_scenarios{2} = tmp_mfreq;

% Scenario 4
tmp_mfreq.crsspctrm = tmp_mfreq.crsspctrm+UU;
mfreq_scenarios{4} = tmp_mfreq;

% Scenario 3
tmp_mfreq = connected_mfreq;
tmp_mfreq.crsspctrm = tmp_mfreq.crsspctrm+UU;
mfreq_scenarios{3} = tmp_mfreq;

%% Power
clf
this_mfreq = mfreq_scenarios{1};
plot(this_mfreq.freq, squeeze(this_mfreq.crsspctrm(1,1,:))); hold on;
plot(this_mfreq.freq, squeeze(this_mfreq.crsspctrm(2,2,:))); 

xlabel('Frequency')
ylabel('Power')
set(gca,'fontsize',20)
legend({'Y1Y1','Y2Y2'});
ylim([0 6]);

%%

% Scenario 3 and 4 are the same
hold on
this_mfreq = mfreq_scenarios{3};
plot(this_mfreq.freq, squeeze(this_mfreq.crsspctrm(1,1,:))); hold on;
plot(this_mfreq.freq, squeeze(this_mfreq.crsspctrm(2,2,:))); 
plot(this_mfreq.freq, this_mfreq.freq*0+UU)

xlabel('Frequency')
ylabel('Power')
set(gca,'fontsize',20)
legend({'Y1Y1','Y2Y2','UU'})
ylim([0 6]);

%%

%chk against simulated time domain data
% empirical_freq.pwrspct=squeeze(mean(abs(empirical_freq_common_input.fourierspctrm).^2,1));
% 
% plot(empirical_freq.freq, cfg.fsample*squeeze(empirical_freq.pwrspct(1,:)) ,'k--'); hold on;
% plot(empirical_freq.freq, cfg.fsample*squeeze(empirical_freq.pwrspct(2,:)),'k--'); 
% plot(empirical_freq.freq, empirical_freq.freq*0+UU)

%% The effect of common input signal on coherence 
% analytic coherence - no common input Scenario 1
this_mfreq = mfreq_scenarios{1};
clf
cfg        = [];
cfg.method = 'coh';
coh = ft_connectivityanalysis(cfg,this_mfreq);
hl_coh_orig=plot(coh.freq,squeeze(coh.cohspctrm(1,2,:)).^2,'b');
ylim([0 0.5]);

%%

%chk
% cfg        = [];     
% cfg.method = 'coh';
% empirical_coh = ft_connectivityanalysis(cfg,empirical_freq);
% hold on
% plot(empirical_coh.freq,squeeze(empirical_coh.cohspctrm(1,2,:)).^2)

% Scenario 2
this_mfreq = mfreq_scenarios{2};
clf
cfg        = [];
cfg.method = 'coh';
coh = ft_connectivityanalysis(cfg,this_mfreq);
hl_coh_orig=plot(coh.freq,squeeze(coh.cohspctrm(1,2,:)).^2,'b');
ylim([0 0.5]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
this_mfreq = mfreq_scenarios{1};
% analytic coherence - common input scenario 3
analyticCohWithInput = @(XiXj,UU,XiXi,XjXj) abs(XiXj+UU).^2./((XiXi+UU).*(XjXj+UU));
XiXi=squeeze(this_mfreq.crsspctrm(1,1,:));
XjXj=squeeze(this_mfreq.crsspctrm(2,2,:));
XiXj=squeeze(this_mfreq.crsspctrm(1,2,:));
clf

analytic_coh_with_input=analyticCohWithInput(XiXj,UU,XiXi,XjXj);
hl_coh_commonInput1=plot(coh.freq,analytic_coh_with_input,'b--');
ylim([0 0.5]);
%%
%chk
% cfg        = [];     
% cfg.method = 'coh';
% empirical_coh = ft_connectivityanalysis(cfg,empirical_freq_common_input);
% hold on
% plot(empirical_coh.freq,squeeze(empirical_coh.cohspctrm(1,2,:)).^2)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analytic coherence - common input scenario 4
this_mfreq = mfreq_scenarios{2};
XiXi=squeeze(this_mfreq.crsspctrm(1,1,:));
XjXj=squeeze(this_mfreq.crsspctrm(2,2,:));
XiXj=squeeze(this_mfreq.crsspctrm(1,2,:));

clf
analytic_coh_with_input=analyticCohWithInput(XiXj,UU,XiXi,XjXj);
hl_coh_commonInput2=plot(coh.freq,analytic_coh_with_input,'b-.');
ylim([0 0.5]);
set([ hl_coh_commonInput2 ],'linewidth',2)
xlabel('Frequency')
ylabel('Coherence')
set(gca,'fontsize',20)

%%

% Ratio identities
NCR1=squeeze(this_mfreq.crsspctrm(1,1,:)/UU);
NCR2=squeeze(this_mfreq.crsspctrm(2,2,:)/UU);

% todo: add equation numbers

% chk, coherence based on ratios
% analyticCohWithInput_Ratio = @(NCR1,NCR2) 1./(NCR1.*NCR2+NCR1+NCR2+1);
% chk_analytic_coh=analyticCohWithInput_Ratio(NCR1,NCR2);
% hold on
% plot(chk_analytic_coh,'r--')

%%%%%%%%%%%%%%%%%%%%%%
%estimated NCR based on coh

plot(coh.freq,NCR1,coh.freq,NCR2);
% estimated CNR
hold on;
est_CNR=1./sqrt(analytic_coh_with_input)-1;
plot(coh.freq,est_CNR,'k--');

xlabel('Frequency')
ylabel('Power')
set(gca,'fontsize',20)
legend({'NCR1','NCR2','NCR'})



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example - common signal introduces inst causality
close all

clf

for scenarioIndex=1:length(mfreq_scenarios)
    
    
    this_mfreq=mfreq_scenarios{scenarioIndex};
    
    % for all but the first scenario we need to set these
    % to zero to make sure they are re-estimated
    if scenarioIndex>1
        this_mfreq=rmfield(this_mfreq, {'transfer','noisecov'});
        
        % if we use Wilson's factorization to estimate the transfer
        % function we have to be careful about the fieldtrip conventions.
        % As per lines 68 - 105 in sfactorization_wilson.m: 
        
        %"
        %check whether the last frequency bin is strictly real-valued.
        % if that's the case, then it is assumed to be the Nyquist frequency
        % and the two-sided spectral density will have an even number of
        % frequency bins. if not, in order to preserve hermitian symmetry,
        % the number of frequency bins needs to be odd."
        
        % "% the input cross-spectral density is assumed to be weighted with a
        % factor of 2 in all non-DC and Nyquist bins, therefore weight the 
        % DC-bin with a factor of 2 to get a correct two-sided
        % representation"
        
        %because fieldtip will assume these, and since we have not applied
        %any weighting, we enforce this now
        this_mfreq.crsspctrm(:,:,1)=this_mfreq.crsspctrm(:,:,1)/2;
        this_mfreq.crsspctrm(:,:,end)=this_mfreq.crsspctrm(:,:,end)/2;

        % if we don't do this we will simply see some ringing due to the
        % applicatios of fft/ifft. 
    end

    
    cfg        = [];
    cfg.method = 'granger';
    GC = ft_connectivityanalysis(cfg,this_mfreq);

    cfg.method = 'total_interdependence';
    total_interdependence = ft_connectivityanalysis(cfg,this_mfreq);

    cfg.method = 'coh';
    coh = ft_connectivityanalysis(cfg,this_mfreq);

    cfg.method = 'instantaneous_causality';
    inst_GC = ft_connectivityanalysis(cfg, this_mfreq);
    
     %chk - relationship between coherence and total interdependence
    % % clf
    % % plot(analytic_total_interdependence.freq,squeeze(analytic_total_interdependence.totispctrm(1,2,:)),'r');hold on;
    % % plot(analytic_coh.freq,-log(1-squeeze(analytic_coh.cohspctrm(1,2,:).^2)),'k--');hold on;
    
    %chk - empirical estimation for scenario 1 and 3 
    if scenarioIndex == 1 || scenarioIndex == 3
        
        if scenarioIndex == 1
            this_mfreq = empirical_freq;
        else
            this_mfreq = empirical_freq_common_input;
        end
        
        cfg        = [];
        cfg.method = 'granger';
        GC_E = ft_connectivityanalysis(cfg,this_mfreq);

        cfg.method = 'total_interdependence';
        total_interdependence_E = ft_connectivityanalysis(cfg,this_mfreq);

        cfg.method = 'coh';
        coh_E = ft_connectivityanalysis(cfg,this_mfreq);

        cfg.method = 'instantaneous_causality';
        inst_GC_E = ft_connectivityanalysis(cfg, this_mfreq);
    end

   

    subplot_counters=scenarioIndex:4:16;

    subplot(4,4,subplot_counters(1))
    h(1)=plot(coh.freq,squeeze(coh.cohspctrm(1,2,:).^2),'b');hold on;
    h(2)=plot(coh.freq,-log(1-squeeze(coh.cohspctrm(1,2,:).^2)),'r','linewidth',3);
    legend('C', '-log(1-C)')
    ylim([0 0.6])

    subplot(4,4,subplot_counters(2))
    h(4)=plot(inst_GC.freq,squeeze(inst_GC.instantspctrm(2,1,:)),'b');
    ylim([0 0.6])

    subplot(4,4,subplot_counters(3))
    h(3)=plot(GC.freq,squeeze(GC.grangerspctrm(1,2,:)+GC.grangerspctrm(2,1,:)),'b');
    ylim([0 0.2])
    
    subplot(4,4,subplot_counters(4))
    precent_inst=100*squeeze(inst_GC.instantspctrm(2,1,:))./squeeze(-log(1-squeeze(coh.cohspctrm(1,2,:).^2)));
    plot(GC.freq,precent_inst,'b','linewidth',2);
    ylim([0 100])
    xlim([0 200])
    
    
    % add empirical check against time domain data (dashed lines)
%     if scenarioIndex == 1 || scenarioIndex == 3
%         subplot(4,4,subplot_counters(1))
%         h(1)=plot(coh.freq,squeeze(coh_E.cohspctrm(1,2,:).^2),'b--');hold on;
%         h(2)=plot(coh.freq,-log(1-squeeze(coh_E.cohspctrm(1,2,:).^2)),'r--','linewidth',3);
%         legend('C', '-log(1-C)')
%         ylim([0 0.6])
%         
%         subplot(4,4,subplot_counters(2))
%         hold on
%         h(4)=plot(inst_GC_E.freq,squeeze(inst_GC_E.instantspctrm(2,1,:)),'b--');
%         ylim([0 0.6])
%         
%         subplot(4,4,subplot_counters(3))
%         hold on
%         h(3)=plot(GC.freq,squeeze(GC_E.grangerspctrm(1,2,:)+GC_E.grangerspctrm(2,1,:)),'b--');
%         ylim([0 0.2])
%         
%         subplot(4,4,subplot_counters(4))
%         hold on
%         precent_inst=100*squeeze(inst_GC_E.instantspctrm(2,1,:))./squeeze(-log(1-squeeze(coh_E.cohspctrm(1,2,:).^2)));
%         plot(GC.freq,precent_inst,'b','linewidth',2);
%     end
    
end


%%


