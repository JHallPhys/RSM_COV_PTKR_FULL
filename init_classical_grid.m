function [q,p,dq,dp,qmesh,pmesh]=init_classical_grid(qL,qR,pL,pR,N)

q=linspace(qL,qR,N); % The initial condition for your q
p=linspace(pL,pR,N); % The range of your p - set to 0 for simplicity can be a range
dq=abs(q(2)-q(1));
dp=abs(p(2)-p(1));
[qmesh,pmesh]=meshgrid(q,p);
end