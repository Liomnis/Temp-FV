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
%% parámetros
Uo = 30.02;
U1 = 6.28;
UL = (G_noct*0.9)/(T_noct-Ta_noct);
% %% interpolar si algún valor medido es cero
% if any(Ta_data==0) %si hay algun valor=0
%     x=1:length(Ta_data);
%     i=find(~(Ta_data==0));
%     Ta_data=interp1(x(i),Ta_data(i),x);
% end
% if any(Tc_data==0) %si hay algun valor=0
%     x=1:length(Tc_data);
%     i=find(~(Tc_data==0));
%     Tc_data=interp1(x(i),Tc_data(i),x);
% end
%% validar solo radiaciones mayores de 400 W/m2
c = 1;
for i=1:length(G_data)
    if G_data(i) >= 1
        G(1,c) = G_data(i);
        Tc(1,c) = Tc_data(i);
        Ta(1,c) = Ta_data(i);
        v(1,c) = v_data(i);
        c = c + 1;
    end
end
%% normalizar
G_normal = G / max(G);
G_simANN = G_normal(1,1:(length(G_normal)*3/4));
G_valANN = G_normal(1,(length(G_normal)*3/4)+1:(length(G_normal)));
Ta_normal = Ta / max(Ta);
Ta_simANN = Ta_normal(1,1:(length(Ta_normal)*3/4));
Ta_valANN = Ta_normal(1,(length(Ta_normal)*3/4)+1:(length(Ta_normal)));
v_normal = v / max(v);
v_simANN = v_normal(1,1:(length(v_normal)*3/4));
v_valANN = v_normal(1,(length(v_normal)*3/4)+1:(length(v_normal)));
Tc_normal = Tc / max(Tc);
% Tc_normal = Tc_normal';
Tc_simANN = Tc_normal(1,1:(length(Tc_normal)*3/4));
Tc_valANN = Tc_normal(1,(length(Tc_normal)*3/4)+1:(length(Tc_normal)));
% 
input_ANNestimado = [Ta_simANN; G_simANN; v_simANN]';
input_ANFISestimado = [Ta_simANN; G_simANN; v_simANN; Tc_simANN]';
output_estimado = Tc_simANN';
%
input_validacion = [Ta_valANN; G_valANN; v_valANN]';
output_validacion = Tc_valANN';
%% estadistica descriptiva
fprintf('\nMinima_Ta= %3.4f    promedio_Ta= %3.4f    Maxima_Ta= %3.4f\n',...
    min(Ta), mean(Ta), max(Ta))
fprintf('\nMinima_Tc= %3.4f    promedio_Tc= %3.4f    Maxima_Tc= %3.4f\n',...
    min(Tc), mean(Tc), max(Tc))
fprintf('\nMinima_G= %3.4f    promedio_G= %3.4f    Maxima_G= %3.4f\n',...
    min(G), mean(G), max(G))
fprintf('\nMinima_v= %3.4f    promedio_v= %3.4f    Maxima_v= %3.4f\n',...
    min(v), mean(v), max(v))
%% cálculo de los modelos
for i=1:length(v)
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
    % 5 Mattei2
    Upv_Mattei2 = 24.1 + 2.9*v(i);
    Tc_Mattei2 = (Upv_Mattei2*Ta(i)+G(i)*(0.81-n_stc*...
        (1-Beta_stc*T_stc)))/(Upv_Mattei2+Beta_stc*n_stc*G(i));
    array_Mattei2(i) = Tc_Mattei2;
    % 6 Tc_Faiman
    Tc_Faiman = Ta(i)+ (G(i)/(Uo + U1*v(i)));
    array_Faiman(i) = Tc_Faiman;
    % 7 Skoplaki2 (2008)
    hw_Sko2 = 8.91 + 2*v(i);
    hw_noct_Sko2 = 8.91 + 2;
    Tc_Sko2 = Ta(i)+(G(i)/G_noct)*(T_noct-Ta_noct)*...
        hw_noct_Sko2/hw_Sko2*(1-n_stc/0.9*1-Beta_stc*T_stc);%
    array_Sko2(i) = Tc_Sko2;
	% 8 Duffie and Beckman (2013)
    Tc_Duffie = Ta(i)+(G(i)/G_noct)*(9.5/(5.7 + 3.8*v(i)))*(T_noct-Ta_noct)*(1-n_stc/0.9);
    array_Duffie(i) = Tc_Duffie;	
	% 9 Risser & Fuentes (1983)
    Tc_Risser = 1.31*Ta(i)+0.0282*G(i)-1.65*v(i)+3.81;
    array_Risser(i) = Tc_Risser;
	% 10 Skoplaki3 (2008)
	Tc_Sko3 = Ta(i)+(0.32/(8.91+2.0*v(i)))*G(i);%
    array_Sko3(i) = Tc_Sko3;
	% 11 Muzathik (2014) 
    Tc_Muzathik = 0.943*Ta(i)+0.0195*G(i)-1.523*v(i)+0.3529 ;
    array_Muzathik(i) = Tc_Muzathik; 
end
%% guardar valore de Tc modeladas en fichero txt
results_Tc_all = [Tc;array_Std;array_Eckstein;array_King;array_Mattei1;...
    array_Mattei2;array_Faiman;array_Sko2;array_Duffie;array_Risser;...
    array_Sko3;array_Muzathik];
fileID = fopen('results_Tc_all.txt','w');
fprintf(fileID,'%8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s\n',...
    'Tc','Std','Eckstein','King','Mattei1','Mattei2','Faiman','Sko2',...
    'Duffie','Risser','Sko3','Muzathik');
fprintf(fileID,'%8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n',...
    results_Tc_all);
fclose(fileID);

results_Tc_data_all = [Tc;array_Std;array_Eckstein;array_King;array_Mattei1;...
    array_Mattei2;array_Faiman;array_Sko2;array_Duffie;array_Risser;...
    array_Sko3;array_Muzathik];
fileID = fopen('results_Tc_data_all.txt','w');
fprintf(fileID,'%8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n',...
    results_Tc_data_all);
fclose(fileID);
%% leer valores de Tc modeladas en fichero txt
medicionesTc = dlmread('results_Tc_data_all.txt');
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
results_Error_all = [modelos;R2_arreglo';RMSE_arreglo';MAE_arreglo'];
fileID = fopen('results_Error_all.txt','w');
fprintf(fileID,'%8s %8s %8s %8s \n', 'modelos','R2','RMSE','MAE');
fprintf(fileID,'%8.0f %8.4f %8.4f %8.4f\n',results_Error_all);
fclose(fileID);
%% graf errores
figure(1)
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
kolor = abs(rand(1,3));
plot(R2_arreglo(11,1), RMSE_arreglo(11,1),'+','Color',(1/255*[233 143 203]))
hold on
xlabel('R^2');
ylabel('RMSE (°C)');
legend('Ross & Smokler (1986)','Eckstein (1990)','King et al.(2004)',...
    'Mattei et al. (2006) v1','Mattei et al. (2006) v2','Faiman (2008)',...
    'Skoplaki et al. (2008)','Duffie & Beckman (2013)',...
    'Risser & Fuentes (1983)','Skoplaki et al. (2009)',...
    'Muzathik (2014)','Location','best')

figure(2)
plot(MAE_arreglo(1,1), RMSE_arreglo(1,1),'s','Color','k','LineWidth',1.5)
hold on
plot(MAE_arreglo(2,1), RMSE_arreglo(2,1),'*','Color',[0.07 0.62 1.00],'LineWidth',1.5)
hold on
plot(MAE_arreglo(3,1), RMSE_arreglo(3,1),'<','Color','b','LineWidth',1.5)
hold on
plot(MAE_arreglo(4,1), RMSE_arreglo(4,1),'>','Color','m','LineWidth',1.5)
hold on
plot(MAE_arreglo(5,1), RMSE_arreglo(5,1),'p','Color','y','LineWidth',1.5)
hold on
plot(MAE_arreglo(6,1), RMSE_arreglo(6,1),'^','Color','g','LineWidth',1.5)
hold on
plot(MAE_arreglo(7,1), RMSE_arreglo(7,1),'o','Color','c','LineWidth',1.5)
hold on
plot(MAE_arreglo(8,1), RMSE_arreglo(8,1),'x','Color',[0.1 0.3250 0.0980],'LineWidth',1.5)
hold on
plot(MAE_arreglo(9,1), RMSE_arreglo(9,1),'d','Color',[1.00 0.41 0.16],'LineWidth',1.5)
hold on
plot(MAE_arreglo(10,1), RMSE_arreglo(10,1),'h','Color',[0 0.5 0],'LineWidth',1.5)
hold on
kolor = abs(rand(1,3));
plot(MAE_arreglo(11,1), RMSE_arreglo(11,1),'+','Color',(1/255*[233 143 203]),'LineWidth',1.5)
hold on
xlabel('MAE (°C)');
ylabel('RMSE (°C)');
legend('Ross & Smokler (1986)','Eckstein (1990)','King et al.(2004)',...
    'Mattei et al. (2006) v1','Mattei et al. (2006) v2','Faiman (2008)',...
    'Skoplaki et al. (2008)','Duffie & Beckman (2013)',...
    'Risser & Fuentes (1983)','Skoplaki et al. (2009)',...
    'Muzathik (2014)','Location','best')

%% graficar errores en barras
figure(3)
x = 1:11;
b = bar(x,round(R2_arreglo, 3, 'significant'));

R2 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(R2,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
axis([0 max(x)+1 0 max(R2_arreglo)+0.5])

%% histogramas
figure(4)
nbins = 10;
Ta_hist = histogram(Ta,nbins);
Ta_hist.FaceColor = 'r';
ylabel('Data observed');
xlabel('Ambient Temp. Ta (°C)');
hold on
%exporta imagen
ax = gca;
exportgraphics(ax,'Ta_Histogram.png','Resolution',60)
%
figure (5)
nbins = 10;
Tc_hist = histogram(Tc,nbins);
Tc_hist.FaceColor = 'k';
ylabel('Data observed');
xlabel('Cell Temp. Tc (°C)');
%exporta imagen
ax = gca;
exportgraphics(ax,'Tc_Histogram.png','Resolution',60)

%% funcion de probabilidad
Ta(find(Ta==0)) = [];
% Extract the unique values occuring in the series
uniqueVals = unique(Ta);
% Get the number of unique values
nbUniqueVals = length(uniqueVals);
% Find the number of occurences of each unique wind speed value
for i=1:nbUniqueVals
    nbOcc = Ta(find(Ta==uniqueVals(i)));
    N(i) = length(nbOcc);
end
% Get the total number of measurements
nbMeas = sum(N);
% To take into account the measurement resolution
% (i.e., a measured wind speed of 2.0 m/s may actually correspond to a
% real wind speed of 2.05 or 1.98 m/s), compute the delta vector which
% contains the difference between two consecutive unique values
delta(1) = uniqueVals(1);
for i=2:(nbUniqueVals)
    delta(i) = uniqueVals(i) - uniqueVals(i-1);
end
% Get the frequency of occurence of each unique value
for i=1:nbUniqueVals
    prob(i) = N(i)/(nbMeas*delta(i));
end
% Get the cumulated frequency
freq = 0;
for i=1:nbUniqueVals
    freq = prob(i)*delta(i) + freq;
    cumFreq(i) = freq;
end

% Plot the distribution Ta
figure(18)
subplot(2,1,1);
plot(uniqueVals,prob,'b');
title('Distribution extracted from the time series');
xlabel('Ambient Temp. Ta (°C)');
ylabel('Probability');
axis([min(uniqueVals) max(uniqueVals)+1 min(prob) max(prob)+0.01])
% Plot the cumulative distribution
subplot(2,1,2);
plot(uniqueVals,cumFreq)
title('Cumulative distribution extracted from the time series');
xlabel('Ambient Temp. Ta (°C)');
ylabel('Cumulative probability');


Tc(find(Tc==0)) = [];
% Extract the unique values occuring in the series
uniqueVals = unique(Tc);
% Get the number of unique values
nbUniqueVals = length(uniqueVals);
% Find the number of occurences of each unique wind speed value
for i=1:nbUniqueVals
    nbOcc = Tc(find(Tc==uniqueVals(i)));
    N(i) = length(nbOcc);
end
% Get the total number of measurements
nbMeas = sum(N);
% To take into account the measurement resolution
% (i.e., a measured wind speed of 2.0 m/s may actually correspond to a
% real wind speed of 2.05 or 1.98 m/s), compute the delta vector which
% contains the difference between two consecutive unique values
delta(1) = uniqueVals(1);
for i=2:(nbUniqueVals)
    delta(i) = uniqueVals(i) - uniqueVals(i-1);
end
% Get the frequency of occurence of each unique value
for i=1:nbUniqueVals
    prob(i) = N(i)/(nbMeas*delta(i));
end
% Get the cumulated frequency
freq = 0;
for i=1:nbUniqueVals
    freq = prob(i)*delta(i) + freq;
    cumFreq(i) = freq;
end

% Plot the distribution Tc
figure(19)
subplot(2,1,1);
plot(uniqueVals,prob,'k');
title('Distribution extracted from the time series');
xlabel('Ambient Temp. Tc (°C)');
ylabel('Probability');
% Plot the cumulative distribution
subplot(2,1,2);
plot(uniqueVals,cumFreq)
title('Cumulative distribution extracted from the time series');
xlabel('Ambient Temp. Tc (°C)');
ylabel('Cumulative probability');


%% cálculo de la eficiencia 
Tc_King = medicionesTc(1:end,4); % Mattei1
Tc_Risser = medicionesTc(1:end,10); % Duffie

nc_best = n_stc*(1 + Beta_stc*(array_Mattei1 - T_stc))*100;
nc_lower = n_stc*(1 + Beta_stc*(array_Duffie - T_stc))*100;

% Graficar ajuste de Eff y Tc
tipoajuste = fittype('poly1');% definir tipo de ajuste
[ajustelinT, bondad] = fit(Tc_King, nc_best', tipoajuste);% se hace el ajuste
p1_King=ajustelinT.p1;% obtener valores de coef del ajuste
p2_King=ajustelinT.p2;% obtener valores de coef del ajuste
R2_King=bondad.rsquare;% obtener valores de coef de correlacion
R2_adj_King=bondad.adjrsquare;% obtener valores de coef de correlacion
RMSE_King=bondad.rmse;
ajusteEff_best = p1_King*Tc_King'+p2_King;

% Duffie
tipoajuste = fittype('poly1');% definir tipo de ajuste
[ajustelinT, bondad] = fit(Tc_Risser, nc_lower', tipoajuste);% se hace el ajuste
p1_Risser=ajustelinT.p1;% obtener valores de coef del ajuste
p2_Risser=ajustelinT.p2;% obtener valores de coef del ajuste
R2_Risser=bondad.rsquare;% obtener valores de coef de correlacion
R2_adj_Risser=bondad.adjrsquare;% obtener valores de coef de correlacion
RMSE_Risser=bondad.rmse;
ajusteEff_lower = p1_Risser*Tc_Risser+p2_Risser;

% grafico ajustes Ta y Tc
figure(6)
plot(Tc_King,ajusteEff_best,'m','LineWidth',1.5)
hold on
plot(Tc_Risser,ajusteEff_lower,'g','LineWidth',1.5)
hold on
xlim([min(Tc_King)-15 max(Tc_King)+15])
ylim([min(ajusteEff_best)-1 max(ajusteEff_lower)+.8])
xlabel('Cell Temp. Tc (°C)')
ylabel('Cell Efficiency  \eta_c(%)')
legend('King et al.(2004)','Risser & Fuentes (1983)','Location','best')
text(10,14,{['\eta_c_ = ',num2str(p2_King),' -',num2str(p1_King),'*T_c']},...
    'HorizontalAlignment','left')
% text(10,13.7,{['\eta_c_ _l_o_w_e_s_t = ',num2str(p2_Risser),' -',num2str(p1_Risser),'*T_c']},...
%     'HorizontalAlignment','left')

%% Graficar ajuste de Ta y Tc
tipoajuste = fittype('poly1');% definir tipo de ajuste
[ajustelinT, bondad] = fit(Ta', array_Std', tipoajuste);% se hace el ajuste
p1=ajustelinT.p1;% obtener valores de coef del ajuste
p2=ajustelinT.p2;% obtener valores de coef del ajuste
R2=bondad.rsquare;% obtener valores de coef de correlacion
R2_adj=bondad.adjrsquare;% obtener valores de coef de correlacion
RMSE=bondad.rmse;
ajusteT_Std = p1*Ta+p2;
% 2
tipoajuste = fittype('poly1');% definir tipo de ajuste
[ajustelinT, bondad] = fit(Ta', array_Eckstein', tipoajuste);% se hace el ajuste
p1=ajustelinT.p1;% obtener valores de coef del ajuste
p2=ajustelinT.p2;% obtener valores de coef del ajuste
R2=bondad.rsquare;% obtener valores de coef de correlacion
R2_adj=bondad.adjrsquare;% obtener valores de coef de correlacion
RMSE=bondad.rmse;
ajusteT_Eckstein = p1*Ta+p2;
% 3
tipoajuste = fittype('poly1');% definir tipo de ajuste
[ajustelinT, bondad] = fit(Ta', array_King', tipoajuste);% se hace el ajuste
p1=ajustelinT.p1;% obtener valores de coef del ajuste
p2=ajustelinT.p2;% obtener valores de coef del ajuste
R2=bondad.rsquare;% obtener valores de coef de correlacion
R2_adj=bondad.adjrsquare;% obtener valores de coef de correlacion
RMSE=bondad.rmse;
ajusteT_King = p1*Ta+p2;
% 4
tipoajuste = fittype('poly1');% definir tipo de ajuste
[ajustelinT, bondad] = fit(Ta', array_Mattei1', tipoajuste);% se hace el ajuste
p1=ajustelinT.p1;% obtener valores de coef del ajuste
p2=ajustelinT.p2;% obtener valores de coef del ajuste
R2=bondad.rsquare;% obtener valores de coef de correlacion
R2_adj=bondad.adjrsquare;% obtener valores de coef de correlacion
RMSE=bondad.rmse;
ajusteT_Mattei1 = p1*Ta+p2;
% 5
tipoajuste = fittype('poly1');% definir tipo de ajuste
[ajustelinT, bondad] = fit(Ta', array_Mattei2', tipoajuste);% se hace el ajuste
p1=ajustelinT.p1;% obtener valores de coef del ajuste
p2=ajustelinT.p2;% obtener valores de coef del ajuste
R2=bondad.rsquare;% obtener valores de coef de correlacion
R2_adj=bondad.adjrsquare;% obtener valores de coef de correlacion
RMSE=bondad.rmse;
ajusteT_Mattei2 = p1*Ta+p2;
% 6
tipoajuste = fittype('poly1');% definir tipo de ajuste
[ajustelinT, bondad] = fit(Ta', array_Faiman', tipoajuste);% se hace el ajuste
p1=ajustelinT.p1;% obtener valores de coef del ajuste
p2=ajustelinT.p2;% obtener valores de coef del ajuste
R2=bondad.rsquare;% obtener valores de coef de correlacion
R2_adj=bondad.adjrsquare;% obtener valores de coef de correlacion
RMSE=bondad.rmse;
ajusteT_Faiman = p1*Ta+p2;
% 7
tipoajuste = fittype('poly1');% definir tipo de ajuste
[ajustelinT, bondad] = fit(Ta', array_Sko2', tipoajuste);% se hace el ajuste
p1=ajustelinT.p1;% obtener valores de coef del ajuste
p2=ajustelinT.p2;% obtener valores de coef del ajuste
R2=bondad.rsquare;% obtener valores de coef de correlacion
R2_adj=bondad.adjrsquare;% obtener valores de coef de correlacion
RMSE=bondad.rmse;
ajusteT_Sko2 = p1*Ta+p2;
% 8
tipoajuste = fittype('poly1');% definir tipo de ajuste
[ajustelinT, bondad] = fit(Ta', array_Duffie', tipoajuste);% se hace el ajuste
p1=ajustelinT.p1;% obtener valores de coef del ajuste
p2=ajustelinT.p2;% obtener valores de coef del ajuste
R2=bondad.rsquare;% obtener valores de coef de correlacion
R2_adj=bondad.adjrsquare;% obtener valores de coef de correlacion
RMSE=bondad.rmse;
ajusteT_Duffie = p1*Ta+p2;
% 9
tipoajuste = fittype('poly1');% definir tipo de ajuste
[ajustelinT, bondad] = fit(Ta', array_Risser', tipoajuste);% se hace el ajuste
p1=ajustelinT.p1;% obtener valores de coef del ajuste
p2=ajustelinT.p2;% obtener valores de coef del ajuste
R2=bondad.rsquare;% obtener valores de coef de correlacion
R2_adj=bondad.adjrsquare;% obtener valores de coef de correlacion
RMSE=bondad.rmse;
ajusteT_Risser = p1*Ta+p2;
% 10 
tipoajuste = fittype('poly1');% definir tipo de ajuste
[ajustelinT, bondad] = fit(Ta', array_Sko3', tipoajuste);% se hace el ajuste
p1=ajustelinT.p1;% obtener valores de coef del ajuste
p2=ajustelinT.p2;% obtener valores de coef del ajuste
R2=bondad.rsquare;% obtener valores de coef de correlacion
R2_adj=bondad.adjrsquare;% obtener valores de coef de correlacion
RMSE=bondad.rmse;
ajusteT_Sko3 = p1*Ta+p2;
% 11
tipoajuste = fittype('poly1');% definir tipo de ajuste
[ajustelinT, bondad] = fit(Ta', array_Muzathik', tipoajuste);% se hace el ajuste
p1=ajustelinT.p1;% obtener valores de coef del ajuste
p2=ajustelinT.p2;% obtener valores de coef del ajuste
R2=bondad.rsquare;% obtener valores de coef de correlacion
R2_adj=bondad.adjrsquare;% obtener valores de coef de correlacion
RMSE=bondad.rmse;
ajusteT_Muzathik = p1*Ta+p2;
% 12
tipoajuste = fittype('poly1');% definir tipo de ajuste
[ajustelinT, bondad] = fit(Ta', Tc, tipoajuste);% se hace el ajuste
p1=ajustelinT.p1;% obtener valores de coef del ajuste
p2=ajustelinT.p2;% obtener valores de coef del ajuste
R2=bondad.rsquare;% obtener valores de coef de correlacion
R2_adj=bondad.adjrsquare;% obtener valores de coef de correlacion
RMSE=bondad.rmse;
ajusteT_Tc = p1*Ta+p2;

% grafico ajustes Ta y Tc
figure(7)
plot(Ta,ajusteT_Tc,'k','LineWidth',2.5)
hold on
plot(Ta,ajusteT_Std,'k')
hold on
plot(Ta,ajusteT_Eckstein,'r')
hold on
plot(Ta,ajusteT_King,'b')
hold on
plot(Ta,ajusteT_Mattei1,'m')
hold on
plot(Ta,ajusteT_Mattei2,'y')
hold on
plot(Ta,ajusteT_Faiman,'g')
hold on
plot(Ta,ajusteT_Sko2,'c')
hold on
plot(Ta,ajusteT_Duffie,'Color',[0.8500 0.3250 0.0980])
hold on
plot(Ta,ajusteT_Risser,'Color',[0.3010 0.7450 0.9330])
hold on
plot(Ta,ajusteT_Sko3,'Color',[0 0.5 0])
hold on
plot(Ta,ajusteT_Muzathik,'Color',(1/255*[233 143 203]))
hold on
xlim([min(Ta) max(Ta)+2])
xlabel('Ambient Temp. Ta (°C)')
ylabel('Cell Temp. Tc (°C)');
legend('Ross & Smokler (1986)','Eckstein (1990)','King et al.(2004)',...
    'Mattei et al. (2006) v1','Mattei et al. (2006) v2','Faiman (2008)',...
    'Skoplaki et al. (2008)','Duffie & Beckman (2013)',...
    'Risser & Fuentes (1983)','Skoplaki et al. (2009)',...
    'Muzathik (2014)','Location','best')
hold all
%% Graficar ajuste de G y Tc
% 1
tipoajuste = fittype('poly1');% definir tipo de ajuste
[ajustelinT, bondad] = fit(G', array_Std', tipoajuste);% se hace el ajuste
p1=ajustelinT.p1;% obtener valores de coef del ajuste
p2=ajustelinT.p2;% obtener valores de coef del ajuste
R2=bondad.rsquare;% obtener valores de coef de correlacion
R2_adj=bondad.adjrsquare;% obtener valores de coef de correlacion
RMSE=bondad.rmse;
ajusteT_Std = p1*G+p2;
% 2
tipoajuste = fittype('poly1');% definir tipo de ajuste
[ajustelinT, bondad] = fit(G', array_Eckstein', tipoajuste);% se hace el ajuste
p1=ajustelinT.p1;% obtener valores de coef del ajuste
p2=ajustelinT.p2;% obtener valores de coef del ajuste
R2=bondad.rsquare;% obtener valores de coef de correlacion
R2_adj=bondad.adjrsquare;% obtener valores de coef de correlacion
RMSE=bondad.rmse;
ajusteT_Eckstein = p1*G+p2;
% 3
tipoajuste = fittype('poly1');% definir tipo de ajuste
[ajustelinT, bondad] = fit(G', array_King', tipoajuste);% se hace el ajuste
p1=ajustelinT.p1;% obtener valores de coef del ajuste
p2=ajustelinT.p2;% obtener valores de coef del ajuste
R2=bondad.rsquare;% obtener valores de coef de correlacion
R2_adj=bondad.adjrsquare;% obtener valores de coef de correlacion
RMSE=bondad.rmse;
ajusteT_King = p1*G+p2;
% 4
tipoajuste = fittype('poly1');% definir tipo de ajuste
[ajustelinT, bondad] = fit(G', array_Mattei1', tipoajuste);% se hace el ajuste
p1=ajustelinT.p1;% obtener valores de coef del ajuste
p2=ajustelinT.p2;% obtener valores de coef del ajuste
R2=bondad.rsquare;% obtener valores de coef de correlacion
R2_adj=bondad.adjrsquare;% obtener valores de coef de correlacion
RMSE=bondad.rmse;
ajusteT_Mattei1 = p1*G+p2;
% 5
tipoajuste = fittype('poly1');% definir tipo de ajuste
[ajustelinT, bondad] = fit(G', array_Mattei2', tipoajuste);% se hace el ajuste
p1=ajustelinT.p1;% obtener valores de coef del ajuste
p2=ajustelinT.p2;% obtener valores de coef del ajuste
R2=bondad.rsquare;% obtener valores de coef de correlacion
R2_adj=bondad.adjrsquare;% obtener valores de coef de correlacion
RMSE=bondad.rmse;
ajusteT_Mattei2 = p1*G+p2;
% 6
tipoajuste = fittype('poly1');% definir tipo de ajuste
[ajustelinT, bondad] = fit(G', array_Faiman', tipoajuste);% se hace el ajuste
p1=ajustelinT.p1;% obtener valores de coef del ajuste
p2=ajustelinT.p2;% obtener valores de coef del ajuste
R2=bondad.rsquare;% obtener valores de coef de correlacion
R2_adj=bondad.adjrsquare;% obtener valores de coef de correlacion
RMSE=bondad.rmse;
ajusteT_Faiman = p1*G+p2;
% 7
tipoajuste = fittype('poly1');% definir tipo de ajuste
[ajustelinT, bondad] = fit(G', array_Sko2', tipoajuste);% se hace el ajuste
p1=ajustelinT.p1;% obtener valores de coef del ajuste
p2=ajustelinT.p2;% obtener valores de coef del ajuste
R2=bondad.rsquare;% obtener valores de coef de correlacion
R2_adj=bondad.adjrsquare;% obtener valores de coef de correlacion
RMSE=bondad.rmse;
ajusteT_Sko2 = p1*G+p2;
% 8
tipoajuste = fittype('poly1');% definir tipo de ajuste
[ajustelinT, bondad] = fit(G', array_Duffie', tipoajuste);% se hace el ajuste
p1=ajustelinT.p1;% obtener valores de coef del ajuste
p2=ajustelinT.p2;% obtener valores de coef del ajuste
R2=bondad.rsquare;% obtener valores de coef de correlacion
R2_adj=bondad.adjrsquare;% obtener valores de coef de correlacion
RMSE=bondad.rmse;
ajusteT_Duffie = p1*G+p2;
% 9
tipoajuste = fittype('poly1');% definir tipo de ajuste
[ajustelinT, bondad] = fit(G', array_Risser', tipoajuste);% se hace el ajuste
p1=ajustelinT.p1;% obtener valores de coef del ajuste
p2=ajustelinT.p2;% obtener valores de coef del ajuste
R2=bondad.rsquare;% obtener valores de coef de correlacion
R2_adj=bondad.adjrsquare;% obtener valores de coef de correlacion
RMSE=bondad.rmse;
ajusteT_Risser = p1*G+p2;
% 10 
tipoajuste = fittype('poly1');% definir tipo de ajuste
[ajustelinT, bondad] = fit(G', array_Sko3', tipoajuste);% se hace el ajuste
p1=ajustelinT.p1;% obtener valores de coef del ajuste
p2=ajustelinT.p2;% obtener valores de coef del ajuste
R2=bondad.rsquare;% obtener valores de coef de correlacion
R2_adj=bondad.adjrsquare;% obtener valores de coef de correlacion
RMSE=bondad.rmse;
ajusteT_Sko3 = p1*G+p2;
% 11
tipoajuste = fittype('poly1');% definir tipo de ajuste
[ajustelinT, bondad] = fit(G', array_Muzathik', tipoajuste);% se hace el ajuste
p1=ajustelinT.p1;% obtener valores de coef del ajuste
p2=ajustelinT.p2;% obtener valores de coef del ajuste
R2=bondad.rsquare;% obtener valores de coef de correlacion
R2_adj=bondad.adjrsquare;% obtener valores de coef de correlacion
RMSE=bondad.rmse;
ajusteT_Muzathik = p1*G+p2;
% 12
tipoajuste = fittype('poly1');% definir tipo de ajuste
[ajustelinT, bondad] = fit(G', Tc, tipoajuste);% se hace el ajuste
p1=ajustelinT.p1;% obtener valores de coef del ajuste
p2=ajustelinT.p2;% obtener valores de coef del ajuste
R2=bondad.rsquare;% obtener valores de coef de correlacion
R2_adj=bondad.adjrsquare;% obtener valores de coef de correlacion
RMSE=bondad.rmse;
ajusteT_Tc = p1*G+p2;
%% grafico ajustes G y Tc
figure(8)
plot(G,ajusteT_Tc,'k','LineWidth',2.5)
hold on
plot(G,ajusteT_Std,'k')
hold on
plot(G,ajusteT_Eckstein,'r')
hold on
plot(G,ajusteT_King,'b')
hold on
plot(G,ajusteT_Mattei1,'m')
hold on
plot(G,ajusteT_Mattei2,'y')
hold on
plot(G,ajusteT_Faiman,'g')
hold on
plot(G,ajusteT_Sko2,'c')
hold on
plot(G,ajusteT_Duffie,'Color',[0.8500 0.3250 0.0980])
hold on
plot(G,ajusteT_Risser,'Color',[0.3010 0.7450 0.9330])
hold on
plot(G,ajusteT_Sko3,'Color',[0 0.5 0])
hold on
plot(G,ajusteT_Muzathik,'Color',(1/255*[233 143 203]))
hold on
xlabel('Irradiance G(W/m^2)')
ylabel('Cell Temp. Tc (°C)');
legend('Ross & Smokler (1986)','Eckstein (1990)','King et al.(2004)',...
    'Mattei et al. (2006) v1','Mattei et al. (2006) v2','Faiman (2008)',...
    'Skoplaki et al. (2008)','Duffie & Beckman (2013)',...
    'Risser & Fuentes (1983)','Skoplaki et al. (2009)',...
    'Muzathik (2014)','Location','best')
hold all