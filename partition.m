function [Norm_out] = partition(Norm_in,eps_fwd_in,eps_bwd_in,set_index_in)

%==========================================================================
% Partition the Norm map 
%==========================================================================

if isequal(set_index_in,'G')
    Norm_in(Norm_in>eps_fwd_in)=NaN;
    Norm_in(~isnan(Norm_in))=0;
    Norm_in(isnan(Norm_in))=1;
    Norm_out=Norm_in;
elseif isequal(set_index_in,'S')
    Norm_in(Norm_in>=eps_bwd_in & Norm_in<=eps_fwd_in)=NaN;
    Norm_in(~isnan(Norm_in))=0;
    Norm_in(isnan(Norm_in))=1;
    Norm_out=Norm_in;
elseif isequal(set_index_in,'L')
    Norm_in(Norm_in<eps_bwd_in)=NaN;
    Norm_in(~isnan(Norm_in))=0;
    Norm_in(isnan(Norm_in))=1;
    Norm_out=Norm_in;
end




end
