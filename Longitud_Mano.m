function [L] = Longitud_Mano(S_meta,End_M,End_P,End_Mi,End_D)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Ld=length(S_meta(:,1));
V_hand=zeros(Ld);V_prox=zeros(Ld);V_med=zeros(Ld);V_dist=zeros(Ld);
for i=1:Ld
    V_hand(i)=norm(S_meta(i,:)-End_M(i,:));
    V_prox(i)=norm(End_M(i,:)-End_P(i,:));
    V_med(i)=norm(End_P(i,:)-End_Mi(i,:));
    V_dist(i)=norm(End_Mi(i,:)-End_D(i,:));
    L=(mean(V_hand(i)+V_prox(i)+V_med(i)+V_dist(i)))/10;    
end
end

