function [SD_av]=get_norm_and_var(qmesh,pmesh,t_start,t_final,Norm_hm_0,k,gamma,N)


SD_av=zeros(N,N); % Initialise the standard deviation array

for itt_time = t_start:t_final

[Norm_hm,qmesh,pmesh] = normbwd(k,gamma,qmesh,pmesh); % Calculate Norm at time t
Norm_hm_0=Norm_hm.*Norm_hm_0; % Update the map
SD = sqrt((Norm_hm-Norm_hm_0).^2);% Calculate standard deviation
% SD = abs(Norm_hm-Norm_hm_0).^2;% Calculate standard deviation
SD_av=SD_av+SD; % Calculate the average


end




end

