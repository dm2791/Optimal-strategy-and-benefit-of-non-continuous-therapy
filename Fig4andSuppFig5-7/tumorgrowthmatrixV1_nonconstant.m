function [ dn ] = tumorgrowthmatrixV1_nonconstant( t,n,g1on,g2,u ,g1off, lambda )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
M=zeros(2,2);
M(1,1)=g1off + (g1on - g1off) * exp(-lambda*t*20);
M(2,2)=g2;
M(2,1)=u;
n=[n(1);n(2)];
dn=M*n;



end

