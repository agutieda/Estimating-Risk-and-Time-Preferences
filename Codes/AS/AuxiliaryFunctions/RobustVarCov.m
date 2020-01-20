function CovMat_theta_hat = RobustVarCov(loglike,x_hat,Y,menuID,hessian,cluster_var,algoList,varargin)
% Function to compute the covariance matrix of the estimator obtained using maximum likelihood

% Number of clusters
N      = length(menuID);
Groups = unique(cluster_var);
M      = length(Groups);
K      = size(hessian,1);

if isempty(varargin)
    estimator = 'robust';
else
    estimator = varargin{1};
end

%% 1) Compute J1

% We get J1 from the Hessian matrix obtained after estimation, multiplied by a minus
J1 = -hessian;

if strcmp(estimator,'robust')
    
    %% 2) Compute J2
    
    % Gradient of the log-likelihood function for each individual
    J2_i = nan(K,M);
    wait_bar_cluster = waitbar(0,'Computing Robust Std. Errors...');
    update_par = 10;
    update_k = round(M/update_par);
    
    tic; flag_time=0;
    for k=1:M
        
        % Find the individuals in a group
        group_idx = (cluster_var == Groups(k));
        nGroup    = sum(group_idx);
        
        % Define likelihood of indiividuals in that group;
        loglike_k = @(x) loglike(x,Y(group_idx),menuID(group_idx),algoList);
        
        % Compute Jacobian of the log-likelihood evaluated at estimated parameter values
        J2_i(:,k) = jacob_fun(loglike_k,x_hat,1e-5).*nGroup;
        
        if mod(k,update_k)==0
            
            if flag_time==0
                time_by_round = toc;
                flag_time = 1;
            end
            
            update_par = update_par-1;
            timeLeft = round(time_by_round*update_par,0);
            
            waitbar(k/M,wait_bar_cluster,['Estimated time left: ' num2str(timeLeft),' seconds']);
            
        end
        
    end
    
    close(wait_bar_cluster);
    
    % Compute J2 as the sum of the Jacobians of each cluster group
    J2 = zeros(K,K);
    for k=1:M
        J2 = J2 + J2_i(:,k)*J2_i(:,k)';
    end
    
    J2 = J2./N;
    
    %% 3) Use J1 and J2 to get the covariance matrix of the estimators
    Asympt_VarMat = J1\J2/J1;
    
    % Divide by N to get finite sample approximation
    CovMat_theta_hat = Asympt_VarMat./N;
    
else
    
    CovMat_theta_hat = inv(J1);
    
end

end

%% Auxiliary Function
function jac_vec = jacob_fun(func,x,step)
% Compute the Jacobian numerically using forward differences

f0 = feval(func,x); n = size(x,1); m = size(f0,1); x0  = x;
jac_vec = zeros(m,n);

for i=1:n
    step2 = step*max(1,x0(i));
    x = x0;
    x(i) = x0(i) + step2;
    jac_vec(1:m,i) = (feval(func,x) - f0)/step2;
end

end

