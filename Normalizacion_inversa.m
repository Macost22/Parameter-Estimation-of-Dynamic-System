function [K11,K12,K22,m] = Normalizacion_inversa(K11_n,K12_n,K22_n,m_n,min7,max7,min810,max810)

Ld=length(K11_n);
K11=zeros(Ld,1);K12=zeros(Ld,1);K22=zeros(Ld,1);m=zeros(Ld,1);
for i=1:Ld
    K11(i)=K11_n(i)*(max810-min810)+min810;
    K12(i)=K12_n(i)*(max810-min810)+min810;
    K22(i)=K22_n(i)*(max810-min810)+min810;
    m(i)=m_n(i)*(max7-min7)+min7;
end
end

