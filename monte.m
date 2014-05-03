function [cdo_data, csq_data] = monte(p_d, para, rr, att, sq_att, runs )
    % We tried to replicate the graphs and data in the paper, so we used a
    % gaussian copula for the simulations, and the number of pools and
    % bonds per pool:
    a_pool = 40; 
    n_pool = 100; 
    
    
    % Create the trenches for the cdo and cdo^2
    n_tran = length(att) - 1;
    sq_tran = find(sq_att(1) == att); 
    absolute = @(notional, payoff, a, b) (b-a) - (max(notional-payoff-a,0) - max(notional-payoff-b,0)); 
    fraction = @(payoff, a, b) payoff / (b-a);
    cdo_pool_payoff = zeros(a_pool, runs); 
   
    % Compute the vector to check if the bond defaulted
    for r = 1:runs     
        uni = coprnd(para(1), a_pool, n_pool );
        default = sum((uni<=p_d)')'; 
        cdo_pool_payoff(:, r) = n_pool-(1-rr)*default;
        if mod(r, 1000) == 0 
        end
    end
    
    % create matrices fior cdo and cdo^2
    cdo_data = zeros(1+n_tran,2);
    csq_data = zeros(1+n_tran,2);
    cdo_tran_default = zeros(a_pool, runs, n_tran);
    cdo_tran_payoff= zeros(a_pool, runs, n_tran);
    cdo_tran_payoff_fraction= zeros(a_pool, runs, n_tran);
    notional = n_pool;
    cdo_pool_loss = notional-cdo_pool_payoff;

    for i=1:n_tran
        cdo_tran_default(:,:,i) = (cdo_pool_loss>(att(i)*notional));
        cdo_tran_payoff(:,:,i) = absolute(notional, cdo_pool_payoff, att(i)*notional, att(i+1)*notional); 
        cdo_tran_payoff_fraction(:,:,i) = fraction(cdo_tran_payoff(:,:,i), att(i)*notional, att(i+1)*notional);
        cdo_data(1+i, 1) = mean(mean(cdo_tran_default(:,:,i)));
        cdo_data(1+i, 2) = mean(mean(cdo_tran_payoff_fraction(:,:,i)));
    end
    
    % compute cd0 - No need to fill the default proba of the global
    % portfolio
    cdo_data(1,1)= NaN; 
    cdo_data(1,2)= mean(mean(cdo_pool_payoff))/notional;
    csq_tran_default= zeros(1, runs, n_tran);
    csq_tran_payoff= zeros(1, runs, n_tran);
    csq_tran_payoff_fraction= zeros(1, runs, n_tran);
    notional = a_pool*n_pool*(sq_att(2) - sq_att(1));
    csq_pool_payoff= sum(cdo_tran_payoff(:,:,sq_tran));
    csq_pool_loss= notional - csq_pool_payoff;
    
    % repeat the same procedure for the cdo^2
    for i=1:n_tran
        csq_tran_default(:,:,i) = (csq_pool_loss>(att(i)*notional));
        csq_tran_payoff(:,:,i) = absolute(notional, csq_pool_payoff, att(i)*notional, att(i+1)*notional); 
        csq_tran_payoff_fraction(:,:,i) = fraction(csq_tran_payoff(:,:,i), att(i)*notional, att(i+1)*notional);
        csq_data(1+i, 1)= mean(mean(csq_tran_default(:,:,i)));
        csq_data(1+i, 2)= mean(mean(csq_tran_payoff_fraction(:,:,i)));
    end
    csq_data(1,1) = NaN; 
    csq_data(1,2)= mean(csq_pool_payoff)/notional;
  
end
