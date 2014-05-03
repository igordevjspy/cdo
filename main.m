%% The Goal here is to reproduce some of the results of the paper The Economics of Structured Finance (Coval, Jurek, Stafford). 
%% The parameters are: 
% - R = Recovery rate
% - rho: pairwise correlation captured through a gaussian copula (othecopulas
% could have been used instead
% - p : default probability
% - number of pools: 40 
% - number of bonds per pool: 100
% - att: attachement and detachement level for the CDO
% - sq_att: attachement and detachement level for the CDO^2

%%
% Case 1: Fixed R, p
rho = 0.2;
p = 0.05;
R = 0.5;
copula = 'gaussian'
att = [0,0.06,0.12,1];
sq_att = [0.06,0.12];

% Plots of the expected payoff as a function of Rho.
Rho = [0:.05:0.99];
len_rho = length(Rho);

junior = zeros(len_rho,4); %wil contain the probability of default, the expected payoff of the cdo and cdo^2 
mezz= zeros(len_rho,4);
senior = zeros(len_rho,4);
glob = zeros(len_rho,4);


for i = 1:len_rho
    r = Rho(i);
    [cdo_data csq_data ] = monte(0.05, r^2, R, [0 .06 .12 1], [.06 .12], 1000);
    
    glob(i,1:2) = cdo_data(1,1:2);
    glob(i,3:4) = csq_data(1,1:2);
    
    junior(i,1:2) = cdo_data(2,1:2);
    junior(i,3:4) = csq_data(2,1:2);
    
    mezz(i,1:2) = cdo_data(3,1);
    mezz(i,3:4) = csq_data(3,1:2);
    
    senior(i,1:2) = cdo_data(4,1:2);
    senior(i,3:4) = csq_data(4,1:2);
    
    
       
end
    
   
plot(Rho,glob(:,2)/glob(5,2)) ; hold all
plot(Rho,junior(:,2)/junior(5,2))
plot(Rho,mezz(:,2)/mezz(5,2))
plot(Rho,senior(:,2)/senior(5,2))
hleg1 = legend('Junior','Mezzanine', 'Senior', 'Global');



%%
% Case 2: Add some noise: 
% case 2.1: p and rho follow beta distributions with
% mean and variance:
% p: mean 0.05, variance 0.03, 
% sqrt(rho): mean 0.2, variance 0.15
% basically for each simulation, we draw two beta distributions, with
% params specified by the users, and then compute the payoffs and take the
% mean
% case 2.2: the recovery rate follow a binomial distribution (or a uniform
% distribution

% params for p
m_p = 0.05;
sigma_p = 0.03;
alpha_p =  ((1 - m_p) / sigma_p^2 - 1 / m_p) * m_p ^ 2;
beta_p = alpha_p * (1 / m_p - 1);

% params for rho: instead of having rho deterministically, you can add some
% noise. 
m_rho = 0.2 ; 
sigma_rho = 0.3; 
alpha_rho =  ((1 - m_rho) / sigma_rho^2 - 1 / m_rho) * m_rho ^ 2;
beta_rho = alpha_rho * (1 / m_rho - 1);


junior = zeros(len_rho,4);
mezz = zeros(len_rho,4);
senior = zeros(len_rho,4);
glob = zeros(len_rho,4);

n_sim = 100;
% the code is only for p varying, but you can easily adapt the for loop and
% add noise. 

% note that this simulation is quite long, so you can change the vector rho
% to be one single value for instancem and then change params accordingly

for i = 1:len_rho
    r = Rho(i);
    junior_monte = zeros(n_sim,4);
    mezz_monte = zeros(n_sim,4);
    senior_monte = zeros(n_sim,4);
    glob_monte = zeros(n_sim,4); 
   for k = 1:n_sim
       
        p_beta = betarnd(alpha_p,beta_p);
        %rho_beta = betarnd(alpha_rho,beta_rho)
        %R = binornd(1,0.5)

        
 
        [ cdo_data csq_data ] = monte(p_beta, r^2, R, [0 .06 .12 1], [.06 .12], 100 );
        
        glob_monte(k,1:2) = cdo_data(1,1:2);
        glob_monte(k,3:4) = csq_data(1,1:2);
       
        junior_monte(k,1:2) = cdo_data(2,1:2);
        junior_monte(k,3:4) = csq_data(2,1:2);
        
        mezz_monte(k,1:2) = cdo_data(3,1:2);
        mezz_monte(k,3:4) = csq_data(3,1:2);
        
        senior_monte(k,1:2) = cdo_data(4,1:2);
        senior_monte(k,3:4) = csq_data(4,1:2);
        
        
   end
   
   junior(i,1:4) = mean(junior_monte(:,1:4));
   mezz(i,1:4) = mean(mezz_monte(:,1:4));
   senior(i,1:4) = mean(senior_monte(:,1:4));
   glob(i,1:4) = mean(glob_monte(:,1:4));


   
   
end

plot(Rho,junior(:,2)/junior(5,2)) ;hold all
plot(Rho,mezz(:,2)/mezz(5,2))
plot(Rho,senior(:,2)/senior(5,2))
plot(Rho,glob(:,2)/glob(5,2))
hleg1 = legend('Junior','Mezzanine', 'Senior', 'Global');







