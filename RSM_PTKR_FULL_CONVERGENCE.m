clear all
close all
clc

%==========================================================================
%==========================================================================
% Paramter decleration
%==========================================================================
%==========================================================================

N=1001;    
k=5; 
gamma=0.001; 

t_start=2; % Start time for norm average
t_final=6; % End time for norm average
t_stable=40; % Time for the stable set to be taken up to

nfq=0.44;
tol=1e-4; % Tolerance to stop calculation
max_itteration=100; % Max number of itterations if tolerance not reached

%==========================================================================
%==========================================================================
% Overestimate the gain set
%==========================================================================
%==========================================================================

% Initialise the system

[q,p,dq,dp,qmesh,pmesh]=init_classical_grid(0,1,-0.5,0.5,N);
Norm_hm_av = zeros(N,N,t_final); 
Norm_hm_0=ones(N,N); % W(q,p,t=0)=1 

% Initial Norm evolution to be discarded

if t_start>1
[Norm_hm_av,Norm_hm_0,qmesh,pmesh]=get_norm(qmesh,pmesh,1,t_start-1,Norm_hm_0,Norm_hm_av,k,gamma);
end

% Now calculate the norm at eatch timestep

[Norm_hm_av,Norm_hm_0,qmesh,pmesh]=get_norm(qmesh,pmesh,t_start,t_final,Norm_hm_0,Norm_hm_av,k,gamma);

% Overestimate the support of the gain set

SSO=zeros(N,N); % Array of the support set

% SSO is the overestiamted support set

for itt_time=t_start:t_final
SSO = SSO+partition(Norm_hm_av(:,:,itt_time)./itt_time,1,1,'G');
end

% Divide by the number of itterations that have contributed

SSO=SSO./(length(t_start:t_final)); % max element possible is 1

%==========================================================================
%==========================================================================
% Now calculate the long-time ordered variance set
%==========================================================================
%==========================================================================

% Reset the system

[q,p,dq,dp,qmesh,pmesh]=init_classical_grid(0,1,-0.5,0.5,N);
Norm_hm_av = zeros(N,N,t_final); 
Norm_hm_0=ones(N,N); % W(q,p,t=0)=1 

% Calculate the standard deviation squared (variance) of the long time set

[SS_var]=get_norm_and_var(qmesh,pmesh,t_start,t_stable,Norm_hm_0,k,gamma,N);

%==========================================================================
%==========================================================================
% Conditioning
%==========================================================================
%==========================================================================

% Okay now what we need to do is vary eps_stable to find the subset of the
% overcomplete set that when removed from the phase space area gives the
% correct number of planck cells.

% Define the lower and upper bounds for the variance of the set

eps_lower=0; % Lower bound
eps_upper=exp(gamma*t_final); % Upper Bound
hbar_half_sqrt=sqrt(1/(4*pi*N)); % Definition of hbar_eff
sigma=N*hbar_half_sqrt; % Rescaling for use with unnormalised imgaussfilt

% Now begin conditioning

for itt=1:max_itteration
    itt
    SS=SSO; % The suport overestimated support set than will be reduced
    SS_eps=SS_var; % Use SS_eps as epsilon conditioned support set
    eps_new=(eps_upper-eps_lower)/2; % Define new epsilon 
    eps_new=eps_new+eps_lower;

    % Partition elements with lowest sigma^2<epsilon_new

    SS_eps(SS_eps<eps_new)=NaN;    % Set elements with lowest to NaN
    SS_eps(~isnan(SS_eps))=0;      % Set elements not equal to NaN=0
    SS_eps(isnan(SS_eps))=1;       % Set element equal to NaN=1

    [i1,i2]=find(SS_eps==1); % Get the coordinates of non-zero elements
 
    % If non-zero element of overestimated set is the same as the epsilon
    % reduced sigma set then remove it

    for j = 1:length(i1)
        if abs(SS(i1(j),i2(j)))>=0  % Common non-zero elemnt
            SS(i1(j),i2(j))=0;      % Set equal to zero
        end
    end

    % Now calculate the classical density 
    
    CD=imgaussfilt(SS,sigma); % Classical Density
    nfc=sum(sum(CD*dq*dp)); % Integrate to get fraction of Planck cells

    % Conditioned density is not convergent and too small
    
    if abs(nfc-nfq)>tol && nfc<nfq % Too much has been removed
       eps_lower=eps_lower;
       eps_upper=eps_new;
    end

    % Conditioned density is not convergent and too small

    if abs(nfc-nfq)>tol && nfc>nfq 
        eps_lower=eps_new;
        eps_upper=eps_upper;
        
    end

    % Conditioned density is convergent

    if abs(nfc-nfq)<tol % Convergence
        'Convergance'
        break
    end

end

figure
imagesc(q,p,CD)
colorbar
colormap(viridis)
set(gca,'YDir','normal')
xlabel('q')
ylabel('p')
caxis([0 1])

figure
imagesc(q,p,ones(N,N)-CD-fliplr(flipud(CD)))
colorbar
colormap(viridis)
set(gca,'YDir','normal')
xlabel('q')
ylabel('p')
caxis([0 1])


