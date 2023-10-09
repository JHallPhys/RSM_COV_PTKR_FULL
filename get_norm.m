function [Norm_hm_av,Norm_hm_0,qmesh,pmesh]=get_norm(qmesh,pmesh,t_start,t_final,Norm_hm_0,Norm_hm_av,k,gamma)

for itt_time = t_start:t_final
itt_time
[Norm_hm,qmesh,pmesh] = normbwd(k,gamma,qmesh,pmesh); % Calculate Norm at time t
Norm_hm_0=Norm_hm.*Norm_hm_0; % Update the map


if itt_time==1
    Norm_hm_av(:,:,itt_time)=Norm_hm_0;
else
    Norm_hm_av(:,:,itt_time) = Norm_hm_av(:,:,itt_time-1) + Norm_hm_0;
end



end

