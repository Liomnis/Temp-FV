%% Liomnis Osorio Laurencio
% scripts para calcular modelos térmicos

% data en el archivo oct2018_abril2019_CIES CON 162 dias
% V (m/s)	Tc (°C)	Ta (°C)	IG (W/m2)	Power DC

%   #     modelos
%   1     Std 
%   2     Eckstein     
%   3     King  
%   4     Mattei1 
%   5     Mattei2  
%   6     Faiman    
%   7     Sko2  
%   8     Duffie  
%   9     Risser   
%   10    Sko3 
%   11    Muzathik

% dia soleado: 43,49,58,70,83,110,119
% dia parcialmente soleado: 2,41,60,62,63,64,100,103,111,129,148
% dia nublado: 61,80,101,116,137,138
%% limpiar espacio de trabajo
clc
clear
%% definir dia (G irregular dia 138) (dia G normal dia )
dia = 49;
%% leer datos medidos
% variables de entrada
Ta_data = dlmread('oct2018_abril2019_CIES_v1.txt', '\t', 'D2..D23320');
G_data = dlmread('oct2018_abril2019_CIES_v1.txt', '\t', 'E2..E23320');
v_data = dlmread('oct2018_abril2019_CIES_v1.txt', '\t', 'B2..B23320');
% variable de salida
Tc_data = dlmread('oct2018_abril2019_CIES_v1.txt', '\t', 'C2..C23320');
%% datos de catálogo STC y NOCT HEE215MA68
G_noct = 800; %W/m2
T_noct = 45; %°C
Ta_noct = 20; %°C
n_stc = 15.03/100;
Beta_stc = -0.34/100;
T_stc = 25;
Pmpp = 250;
G_stc = 1000;
%% parámetros
Uo = 30.02;
U1 = 6.28;
UL = (G_noct*0.9)/(T_noct-Ta_noct);
%% interpolar si algún valor medido es cero
if any(Ta_data==0) %si hay algun valor=0
    x=1:length(Ta_data);
    i=find(~(Ta_data==0));
    Ta_data=interp1(x(i),Ta_data(i),x);
end
if any(Tc_data==0) %si hay algun valor=0
    x=1:length(Tc_data);
    i=find(~(Tc_data==0));
    Tc_data=interp1(x(i),Tc_data(i),x);
end

%% definir rango de datos del dia
% entrada
time_dia=linspace(0,1,144);% definir 1 dia en tiempo
Ta_dia = Ta_data(144*dia:144*dia+143);
G_dia = G_data(144*dia:144*dia+143);
v_dia = v_data(144*dia:144*dia+143);
% salida
Tc_dia = Tc_data(144*dia:144*dia+143);
%% validar solo radiaciones mayores de x W/m2
G = [];
for i=1:length(G_dia)
    if G_dia(i) >= 1
        G(1,length(G)+1) = G_dia(i);
        Tc(length(G),1) = Tc_dia(i);
        Ta(1,length(G)) = Ta_dia(i);
        v(1,length(G)) = v_dia(i);
        time(1,length(G)) = time_dia(i);
    end
end
%% normalizar
G_normal = G / max(G);
Ta_normal = Ta / max(Ta);
v_normal = v / max(v);
Tc_normal = Tc / max(Tc);
% Tc_normal = Tc_normal';
% 
input_estimado = [Ta_normal; G_normal; v_normal]';

%% ANN
ANNTc = ANNTc_CIES(input_estimado)*67;

%% estadistica descriptiva
fprintf('\nMinima_Ta= %3.4f    promedio_Ta= %3.4f    Maxima_Ta= %3.4f\n',...
    min(Ta), mean(Ta), max(Ta))
fprintf('\nMinima_Tc= %3.4f    promedio_Tc= %3.4f    Maxima_Tc= %3.4f\n',...
    min(Tc), mean(Tc), max(Tc))
fprintf('\nMinima_G= %3.4f    promedio_G= %3.4f    Maxima_G= %3.4f\n',...
    min(G), mean(G), max(G))
fprintf('\nMinima_v= %3.4f    promedio_v= %3.4f    Maxima_v= %3.4f\n',...
    min(v), mean(v), max(v))
%% graficar en el mismo plot G y Ta
figure(1)
colororder({'k','b'})
yyaxis left
plot(time,Ta, time,v)
ylabel('Ambient Temp. Ta (°C)');
axis([0.25 0.85 min(Ta) max(Ta)+3])
yyaxis right
plot(time,G)
datetick('x',15)
xlabel('Time (Hours)')
ylabel('Irradiance G(W/m^2)')
axis([min(time-0.01) max(time+0.01) min(G) max(G)+100])
%exporta imagen
ax = gca;
exportgraphics(ax,'figure1.png','Resolution',50)
%%
figure(2)
colororder({'k','r'})
yyaxis left
plot(time,Ta)
xlabel('Time (Hours)')
ylabel('Ambient Temp. Ta (°C)')
axis([0.1 0.55 min(Ta) max(Ta)+2])
yyaxis right
plot(time,v)
datetick('x',15)
ylabel('v_w (m/s)');
% axis([0.25 0.85 min(v) max(v)+1])
%exporta imagen
ax = gca;
exportgraphics(ax,'figure2.png','Resolution',50)
%% cálculo de los modelos
for i=1:length(Tc)
    % 1 Ross & Smokler - (Ross and Smokler, 1986) 
    Tc_Std = Ta(i)+(G(i)/G_noct)*(T_noct-Ta_noct);%
    array_Std(i) = Tc_Std;
	% 2 Eckstein (1990)
	Tc_Eckstein = Ta(i)+(G(i)*0.9/UL)*(1-n_stc/0.9);
    array_Eckstein(i) = Tc_Eckstein;  
	% 3 King et al. (2004) 	
	Tc_King = Ta(i)+G(i)*exp(-3.473-0.0594*v(i));
    array_King(i) = Tc_King;
	% 4 Mattei1
    Upv_Mattei1 = 26.6 + 2.3*v(i);
    Tc_Mattei1 = (Upv_Mattei1*Ta(i)+G(i)*(0.81-n_stc*...
        (1-Beta_stc*T_stc)))/(Upv_Mattei1+Beta_stc*n_stc*G(i));
    array_Mattei1(i) = Tc_Mattei1;
    % 4 Mattei2
    Upv_Mattei2 = 24.1 + 2.9*v(i);
    Tc_Mattei2 = (Upv_Mattei2*Ta(i)+G(i)*(0.81-n_stc*...
        (1-Beta_stc*T_stc)))/(Upv_Mattei2+Beta_stc*n_stc*G(i));
    array_Mattei2(i) = Tc_Mattei2;
    % 5 Tc_Faiman
    Tc_Faiman = Ta(i)+ (G(i)/(Uo + U1*v(i)));
    array_Faiman(i) = Tc_Faiman;
    % 6 Skoplaki2 (2008)
    hw_Sko2 = 5.7 + 2.8*v(i);
    hw_noct_Sko2 = 5.7 + 2.8;
    Tc_Sko2 = Ta(i)+(G(i)/G_noct)*(T_noct-Ta_noct)*...
        hw_noct_Sko2/hw_Sko2*(1-n_stc/0.9*1-Beta_stc*T_stc);%
    array_Sko2(i) = Tc_Sko2;
	% 7 Duffie and Beckman (2013)
    Tc_Duffie = Ta(i)+(G(i)/G_noct)*(9.5/(5.7 + 3.8*v(i)))*(T_noct-Ta_noct)*(1-n_stc/0.9);
    array_Duffie(i) = Tc_Duffie;	
	% 8 Risser & Fuentes (1983)
    Tc_Risser = 1.31*Ta(i)+0.0282*G(i)-1.65*v(i)+3.81;
    array_Risser(i) = Tc_Risser;
	% 9 Skoplaki3 (2008)
	Tc_Sko3 = Ta(i)+(0.32/(5.7+2.8*v(i)))*G(i);%
    array_Sko3(i) = Tc_Sko3;
	% 10 Muzathik (2014) 
    Tc_Muzathik = 0.943*Ta(i)+0.0195*G(i)-1.523*v(i)+0.3529 ;
    array_Muzathik(i) = Tc_Muzathik; 
end
%% guardar valores de Tc modeladas en fichero txt
results_Tc_dia = [Tc';array_Std;array_Eckstein;array_King;array_Mattei1;...
    array_Mattei2;array_Faiman;array_Sko2;array_Duffie;array_Risser;...
    array_Sko3;array_Muzathik];
fileID = fopen('results_Tc_dia.txt','w');
fprintf(fileID,'%8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s\n',...
    'Tc','Std','Eckstein','King','Mattei1','Mattei2','Faiman','Sko2',...
    'Duffie','Risser','Sko3','Muzathik');
fprintf(fileID,'%8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n',...
    results_Tc_dia);
fclose(fileID);

results_Tc_data_dia = [Tc';array_Std;array_Eckstein;array_King;array_Mattei1;...
    array_Mattei2;array_Faiman;array_Sko2;array_Duffie;array_Risser;...
    array_Sko3;array_Muzathik];
fileID = fopen('results_Tc_data_dia.txt','w');
fprintf(fileID,'%8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n',...
    results_Tc_data_dia);
fclose(fileID);

%% graficar en el mismo plot G y Tc
figure(3)
plot(time,array_Std,'k')
hold on
plot(time,array_Eckstein,'Color',[0.07 0.62 1.00])
hold on
plot(time,array_King,'b')
hold on
plot(time,array_Mattei1,'m')
hold on
plot(time,array_Mattei2,'y')
hold on
plot(time,array_Faiman,'g')
hold on
plot(time,array_Sko2,'c')
hold on
plot(time,array_Duffie,'Color',[0.1 0.3250 0.0980])
hold on
plot(time,array_Risser,'Color',[1.00 0.41 0.16])
hold on
plot(time,array_Sko3,'Color',[0 0.5 0])
hold on
plot(time,array_Muzathik,'Color',(1/255*[233 143 203]))
hold on
plot(time,Tc,'k-.')
hold on
plot(time,Ta,'b-.')
hold on
ylabel('Cell Temp. Tc (°C)')
yyaxis right
plot(time,G)
datetick('x',15)
xlabel('Time (Hours)')
ylabel('Irradiance G(W/m^2)')
axis([min(time-0.01) max(time+0.01) min(G) max(G)+100])
legend({'Ross & Smokler (1986)','Eckstein (1990)','King et al.(2004)',...
    'Mattei et al. (2006) v1','Mattei et al. (2006) v2','Faiman (2008)',...
    'Skoplaki et al. (2008)','Duffie & Beckman (2013)',...
    'Risser & Fuentes (1983)','Skoplaki et al. (2009)',...
    'Muzathik (2014)','Measured Tc','Measured Ta','Measured G'},'Location','best')

%exporta imagen
ax = gca;
exportgraphics(ax,'figure3.png','Resolution',50)
%% graficar en el mismo plot G y Tc
figure(4)
colororder({'m','g','k','r','b'})
plot(time,array_King,'LineWidth',1.5)
hold on
plot(time,array_Risser,'LineWidth',1.5)
hold on
plot(time,Tc','.-.')
hold on
plot(time,Ta,'.-.')
hold on
ylabel('Cell Temp. Tc (°C)')
yyaxis right
plot(time,G,'.-.')
datetick('x',15)
xlabel('Time (Hours)')
ylabel('Irradiance G(W/m^2)')
axis([min(time-0.01) max(time+0.01) min(G) max(G)+100])
legend({'Tc_ King et al.(2004)','Tc_ Risser & Fuentes (1983)','Tc_M_e_a_s_u_r_e_d',...
    'Ta_M_e_a_s_u_r_e_d','G_M_e_a_s_u_r_e_d'},'Location','best')
%exporta imagen
ax = gca;
exportgraphics(ax,'figure4.png','Resolution',50)
%% guardar valores de Tc modeladas en fichero txt
medicionesTc = dlmread('results_Tc_data_dia.txt');
Tc = medicionesTc(1:end,1);
%% calcular errores
for i=2:12
    Tc_v = medicionesTc(1:end,i);
    tipoajuste = fittype('poly1');% definir tipo de ajuste
    [ajustelinT, bondad] = fit(Tc_v, Tc, tipoajuste);% se hace el ajuste
    p1=ajustelinT.p1;% obtener valores de coef del ajuste
    p2=ajustelinT.p2;% obtener valores de coef del ajuste
    R2=bondad.rsquare;% obtener valores de coef de correlacion
    R2_adj=bondad.adjrsquare;% obtener valores de coef de correlacion
    RMSE=bondad.rmse;
    
    R2_arreglo(i-1,1) = R2;
    RMSE_arreglo(i-1,1) = RMSE;
    
    MAE_arreglo(i-1,1) = mean(abs(Tc_v-Tc));
end
modelos=linspace(1,length(R2_arreglo),length(R2_arreglo));
results_Error_dia = [modelos;R2_arreglo';RMSE_arreglo';MAE_arreglo'];
fileID = fopen('results_Error_dia.txt','w');
fprintf(fileID,'%8s %8s %8s %8s \n', 'modelos','R2','RMSE','MAE');
fprintf(fileID,'%8.0f %8.4f %8.4f %8.4f\n',results_Error_dia);
fclose(fileID);
%% graf errores
figure(5)
plot(R2_arreglo(1,1), RMSE_arreglo(1,1),'s','Color','k')
hold on
plot(R2_arreglo(2,1), RMSE_arreglo(2,1),'*','Color',[0.07 0.62 1.00])
hold on
plot(R2_arreglo(3,1), RMSE_arreglo(3,1),'<','Color','b')
hold on
plot(R2_arreglo(4,1), RMSE_arreglo(4,1),'>','Color','m')
hold on
plot(R2_arreglo(5,1), RMSE_arreglo(5,1),'p','Color','y')
hold on
plot(R2_arreglo(6,1), RMSE_arreglo(6,1),'^','Color','g')
hold on
plot(R2_arreglo(7,1), RMSE_arreglo(7,1),'o','Color','c')
hold on
plot(R2_arreglo(8,1), RMSE_arreglo(8,1),'x','Color',[0.1 0.3250 0.0980])
hold on
plot(R2_arreglo(9,1), RMSE_arreglo(9,1),'d','Color',[1.00 0.41 0.16])
hold on
plot(R2_arreglo(10,1), RMSE_arreglo(10,1),'h','Color',[0 0.5 0])
hold on
plot(R2_arreglo(11,1), RMSE_arreglo(11,1),'+','Color',(1/255*[233 143 203]))
hold on
xlabel('R^2');
ylabel('RMSE (°C)');
legend('Ross & Smokler','Eckstein','King','Mattei1','Mattei2','Faiman','Skoplaki1',...
    'Duffie & Beckman' ,'Risser & Fuentes ','Skoplaki2',...
    'Muzathik','Location','best')
%exporta imagen
ax = gca;
exportgraphics(ax,'Fit_R2_RMSE_Scatters.png','Resolution',60)

%% cálculo de la eficiencia 
Tc_Risser = medicionesTc(1:end,10); % 
Tc_King = medicionesTc(1:end,4); % 
Tc = Tc';
nc_Risser = n_stc*(1 + Beta_stc*(array_Risser - T_stc))*100;
nc_King = n_stc*(1 + Beta_stc*(array_King - T_stc))*100;
nc_med = n_stc*(1 + Beta_stc*(Tc - T_stc))*100;
%
figure(6)
plot(time,nc_King,'m',time,nc_Risser,'g',time,nc_med,'k.-.','LineWidth',1.5)
hold on
xlabel('Time (Hours)')
ylabel('Cell Efficiency  n_c (%)')
datetick('x',15)
legend('\eta_c King et al.(2004)','\eta_c Risser & Fuentes (1983)',...
    '\eta_c Measured','Location','best')
%exporta imagen
ax = gca;
exportgraphics(ax,'Eff_Line.png','Resolution',60)


% %% Graficar ajuste de Eff y Tc
% tipoajuste = fittype('poly1');% definir tipo de ajuste
% [ajustelinT, bondad] = fit(Tc_King, nc_King', tipoajuste);% se hace el ajuste
% p1=ajustelinT.p1;% obtener valores de coef del ajuste
% p2=ajustelinT.p2;% obtener valores de coef del ajuste
% R2=bondad.rsquare;% obtener valores de coef de correlacion
% R2_adj=bondad.adjrsquare;% obtener valores de coef de correlacion
% RMSE=bondad.rmse;
% ajusteEff_best = p1*Tc_King'+p2;
% 
% % Duffie
% tipoajuste = fittype('poly1');% definir tipo de ajuste
% [ajustelinT, bondad] = fit(Tc_Risser, nc_Risser', tipoajuste);% se hace el ajuste
% p1=ajustelinT.p1;% obtener valores de coef del ajuste
% p2=ajustelinT.p2;% obtener valores de coef del ajuste
% R2=bondad.rsquare;% obtener valores de coef de correlacion
% R2_adj=bondad.adjrsquare;% obtener valores de coef de correlacion
% RMSE=bondad.rmse;
% ajusteEff_lowest = p1*Tc_Risser+p2;
% 
% % grafico ajustes Tc y n_stc
% figure(7)
% plot(Tc_King,ajusteEff_best,'m','LineWidth',1.5)
% hold on
% plot(Tc_Risser,ajusteEff_lowest,'g','LineWidth',1.5)
% hold on
% xlim([min(Tc_King)-4 max(Tc_King)+4])
% ylim([min(ajusteEff_best)-.2 max(ajusteEff_best)+.2])
% xlabel('Tc (°C)')
% ylabel('n_c (%)')
% legend('King','Risser & Fuentes','Location','best')
% text(min(Tc_King)-1,max(ajusteEff_best)-.2,{['n_c = ',num2str(p2),...
%     ' ',num2str(p1),'*T_c']},...
%     'HorizontalAlignment','left')
% 
% %exporta imagen
% ax = gca;
% exportgraphics(ax,'figure7.png','Resolution',50)
