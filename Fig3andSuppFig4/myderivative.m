function dXdt = myderivative(N,M)

S = N(1);
R = N(2);
pgrowth1=M(1,1);
pgrowth2=M(2,2);
mu=M(2,1);

% dCdt = p.growth*C-p.death*C;
% dSdt = p.growth*S-p.death*S;
dSdt = pgrowth1*S;
dRdt = pgrowth2*R+mu*S;

dXdt = [dSdt dRdt];