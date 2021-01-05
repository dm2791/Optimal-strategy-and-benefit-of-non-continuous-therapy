function dXdt = myderivative4celltypes(N,M)

DS = N(1);
RA = N(2);
RB = N(3);
DR = N(4);
pgrowth1=M(1,1);
pgrowth2=M(2,2);
pgrowth3=M(3,3);
pgrowth4=M(4,4);
muDS_A=M(2,1);
muDS_B=M(3,1);
muA_DR=M(4,2);
muB_DR=M(4,3);


dDSdt = pgrowth1*DS;
dRAdt = pgrowth2*RA+muDS_A*DS;
dRBdt = pgrowth3*RB+muDS_B*DS;
dDRdt = pgrowth4*DR+muA_DR*RA+muB_DR*RB;

dXdt = [dDSdt dRAdt dRBdt dDRdt];