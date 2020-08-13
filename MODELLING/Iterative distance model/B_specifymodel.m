function model = C_specifymodel(modelidx)

% Model 1: One iteration,  with k_sig_scale and b_sig_scale
         % starting point is cursor spawn location

switch modelidx
    case {1}
        f_NLL = @get_NLL;
        f_simulate = @func_iter_avg_gaussian_oneiter_bk;
        pars = {'(1) k_sig_scale: sigma distance-scaling factor',...
                '(2) b_sig: sigma baseline',...
                '(3) alpha: mapping between objective and subjective reward',...
                '(4) beta: inverse temperature for confidence resp, higher = lower noise'};
        npars = length(pars);
        lb = [log(0),... %k_sig_scale
              log(0),... %b_sig
              log(0),... %alpha
              log(0)];%,... %beta
        ub = [inf,... %k_sig_scale
              inf,... %b_sig
              inf,... %alpha
              inf];%,... %beta
        plb = [log(0.01),... %k_sig_scale
              log(0.01),... %b_sig
              log(0.01),... %alpha
              log(1)];%,... %beta
        pub = [log(1),... %k_sig_scale
              log(1),... %b_sig
              log(0.1),... %alpha
              log(20)];%,... %beta    
end

model.npars           = npars;
model.lb              = lb;
model.ub              = ub;
model.plb             = plb;
model.pub             = pub;

model.f_NLL           = f_NLL;
model.f_simulate      = f_simulate;