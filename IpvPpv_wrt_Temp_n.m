clear all
clc
Va=0:.01:22;
Suns=1;
TaC=10:5:30;
lva=length(Va);
lT=length(TaC);

Ipv=zeros(size(Va));
for s=1:1:lT    
for i=1:1:lva
k=1.38e-23;
q=1.6e-19;
A=1.2;
Vg=1.12;
Ns=36;
T1=273+25;
Voc_T1=21.06/Ns;
Isc_T1=3.80;
T2=273+75;
Voc_T2=17.05/Ns;
Isc_T2=3.92;
TarK=273+TaC(s);
Tref=273+25;
Iph_T1=Isc_T1*Suns;
a=(Isc_T2-Isc_T1)/Isc_T1*1/(T2-T1);
Iph=Iph_T1*(1+a*(TarK-T1));
Vt_T1=k*T1/q;
Ir_T1=Isc_T1/(exp(Voc_T1/(A*Vt_T1))-1);
Ir_T2=Isc_T2/(exp(Voc_T2/(A*Vt_T1))-1);
b=Vg*q/(A*k);
Ir=Ir_T1*(TarK/T1).^(3/A).*exp(-b.*(1./TarK-1/T1));
X2v=Ir_T1/(A*Vt_T1)*exp(Voc_T1/(A*Vt_T1));
dVdI_Voc=-1.15/Ns/2;
Rs=-dVdI_Voc-1/X2v;
Vt_Ta=A*k*TarK/q;
Vc=Va(i)/Ns;
Ia=zeros(size(Vc));
for j=1:1:100
    Ia=Ia-(Iph-Ia-Ir*(exp((Vc+Ia*Rs)/Vt_Ta)-1))./(-1-Ir*(exp((Vc+Ia*Rs)/Vt_Ta)-1).*Rs/Vt_Ta);
end
Ipv(s,i)=Ia;
Ppv(s,i)=Va(i)*Ia;
end
end
axes1 = axes('Parent',figure,'OuterPosition',[0 0.5 1 0.5]);
xlim(axes1,[0 23]);
ylim(axes1,[0 5]);
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');
title('I-V charateristics');
xlabel('V_p_v (V)');
ylabel('I_p_v (A)');
plot1 = plot(Va(1,:),Ipv(:,:),'Parent',axes1,'LineWidth',1.5);
set(plot1(1),'DisplayName','10C T');
set(plot1(2),'DisplayName','15C T');
set(plot1(3),'DisplayName','20C T');
set(plot1(4),'DisplayName','25C T');
set(plot1(5),'DisplayName','30C T');
axes2 = axes('OuterPosition',[0 0 1 0.5]);
xlim(axes2,[0 23]);
ylim(axes2,[0 70]);
box(axes2,'on');
grid(axes2,'on');
hold(axes2,'all');
title('P-V charateristics');
xlabel('V_p_v (V)');
ylabel('P_p_v (W)');
plot2 = plot(Va(1,:),Ppv(:,:),'Parent',axes2,'LineWidth',1.5);
set(plot2(1),'DisplayName','10C T');
set(plot2(2),'DisplayName','15C T');
set(plot2(3),'DisplayName','20C T');
set(plot2(4),'DisplayName','25C T');
set(plot2(5),'DisplayName','30C T');
legend1 = legend(axes2,'show');

legend2 = legend(axes1,'show');
