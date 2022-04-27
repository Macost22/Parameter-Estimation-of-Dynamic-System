%% Datos del Leap Motion

%Input1 = readtable('rotacionnn.csv','Delimiter',',','ReadVariableNames', true);
Input1 = readtable('feaa2.csv','Delimiter',',','ReadVariableNames', true);
%Datos para entrenar

Datos = table2array(Input1);

%% Grafica de la mano llamando la función plot_hand
%función para graficar la mano, las variables de entrada son los datos en
%formato table y el frame en el que se desea iniciar, el paso y el frame en
%el que se desea finalizar -->plot_hand(Input,Frame_start,Step,Frame_end)
%  figure(1)
%  plot_hand(Input1,1,1,2);


%% Definición de variables 
%Longitud datos
Ld=length(Datos(:,1));
%Manitud eslabón
L_L1=zeros(Ld,1); L_L2=zeros(Ld,1);
% Tiempos
Tiempo=zeros(Ld,1); DeltaT=zeros(Ld,1);
% Matrices de rotación
R_Arm=zeros(3,3,Ld); R_Wrist=zeros(3,3,Ld); R_Hand=zeros(3,3,Ld); R_M=zeros(3,3,Ld);R_N=zeros(3,3,Ld);
R_Arm_f=zeros(3,3,Ld); R_Wrist_f=zeros(3,3,Ld); R_Hand_f=zeros(3,3,Ld);

%% Vector para graficar el tiempo
% Cálculo del tiempo
t=table2array(Input1(:,{'Tiempo'}));
for i=1:length(t)
   %se pasa de microsegundos a segundos 
  Tiempo(i)= (t(i)-t(1))/1000000; 
end
Delta_t=mean(diff(Tiempo));
%% Se definen las posiciones x,y,z de los vectores base
%Articulación falange distal, dedo medio, posiciones finales x,y,z.

Arm_xBasis=table2array(Input1(:,{'arm_basis_x_x','arm_basis_x_y','arm_basis_x_z'}));
Arm_yBasis=table2array(Input1(:,{'arm_basis_y_x','arm_basis_y_y','arm_basis_y_z'}));
Arm_zBasis=table2array(Input1(:,{'arm_basis_z_x','arm_basis_z_y','arm_basis_z_z'}));

Hand_xBasis=table2array(Input1(:,{'Hand_basis_x_x','Hand_basis_x_y','Hand_basis_x_z'}));
Hand_yBasis=table2array(Input1(:,{'Hand_basis_y_x','Hand_basis_y_y','Hand_basis_y_z'}));
Hand_zBasis=table2array(Input1(:,{'Hand_basis_z_x','Hand_basis_z_y','Hand_basis_z_z'}));

%Posición x,y,z -> del codo (E), muñeca(W),posición final de la falange distal
%del dedo medio (M_D), posición final del metacarpiano del dedo medio (M_M)

W=table2array(Input1(:,{'wrist_position_x','wrist_position_y','wrist_position_z'}));
M_end=table2array(Input1(:,{'Middle_fin_meta_end_x','Middle_fin_meta_end_y','Middle_fin_meta_end_z'}));
P_end=table2array(Input1(:,{'Middle_fin_prox_end_x','Middle_fin_prox_end_y','Middle_fin_prox_end_z'}));
Mi_end=table2array(Input1(:,{'Middle_fin_inter_end_x','Middle_fin_inter_end_y','Middle_fin_inter_end_z'}));
D_end=table2array(Input1(:,{'Middle_fin_dist_end_x','Middle_fin_dist_end_y','Middle_fin_dist_end_z'}));

Hand_0=table2array(Input1(:,{'hand_position_x','hand_position_y','hand_position_z'}));

%% Análisis de frecuencia
% %Frecuencia de muestreo
 Fs=1/Delta_t;
% T=1/Fs;
% L=length(Tiempo);
% t= (0:L-1).*T;
% Y=fft(Hand_xBasis(:,2));
% P1=2.*(abs(Y(1:L/2)/L));
% f=Fs.*(0:(L/2)-1)./L;
% plot(f,P1)


%% Filtro Butterworth paso bajo
%Orden del filtro
n=3;
%Frecuencia de corte
Fc=0.25; %para los datos feaa2
%Fc=0.4; %frecuencia de corte para rotacionnn
%Frecuencia de corte del filtro normalizado
wn= Fc/(Fs/2);
%Filtro donde se saca los coeficientes de la función de transferencia
[num,den] = butter(n,wn,'Low');

%%

Hand_xBasis_f=zeros(Ld,3);Hand_yBasis_f=zeros(Ld,3);Hand_zBasis_f=zeros(Ld,3);
Arm_xBasis_f=zeros(Ld,3);Arm_yBasis_f=zeros(Ld,3);Arm_zBasis_f=zeros(Ld,3);Hand_0_f=zeros(Ld,3);

for k=1:3
Hand_xBasis_f(:,k)=filtfilt(num,den,Hand_xBasis(:,k));
Hand_yBasis_f(:,k)=filtfilt(num,den,Hand_yBasis(:,k));
Hand_zBasis_f(:,k)=filtfilt(num,den,Hand_zBasis(:,k));

Arm_xBasis_f(:,k)=filtfilt(num,den,Arm_xBasis(:,k));
Arm_yBasis_f(:,k)=filtfilt(num,den,Arm_yBasis(:,k));
Arm_zBasis_f(:,k)=filtfilt(num,den,Arm_zBasis(:,k));
Hand_0_f(:,k)=filtfilt(num,den,Hand_0(:,k));
end
%Grafica del angulo
figure(2)
plot(Tiempo,Arm_zBasis_f(:,1),'r'),hold on,plot(Tiempo,Arm_zBasis(:,1),'c'),
hold off,legend('Velocidad x','Velocidad y'),xlabel('Tiempo [s]'),ylabel('Velocidad [m/s]'),title('Gráfica de la velocidad del centro de la mano')
grid on;

%% Cinemática 3D
% Calculo de los angulos Cardan-Euler secuencia XYZ
for i=1:Ld
   %MATRIZ DE ROTACIÓN DEL BRAZO RESPECTO AL SISTEMA DE REFERENCIAS GLOBAL
   %sin filtro
     R_Arm(:,:,i)=[Arm_xBasis(i,:);Arm_yBasis(i,:);Arm_zBasis(i,:)];
   %con filtr
     R_Arm_f(:,:,i)=[Arm_xBasis_f(i,:);Arm_yBasis_f(i,:);Arm_zBasis_f(i,:)];
     
   %MATRIZ DE ROTACIÓN DE LA MANO RESPECTO AL SISTEMA DE REFERENCIA GLOBAL
   %sin filtro
     R_Hand(:,:,i)=[Hand_xBasis(i,:);Hand_yBasis(i,:);Hand_zBasis(i,:)];
   %con filtro
     R_Hand_f(:,:,i)=[Hand_xBasis_f(i,:);Hand_yBasis_f(i,:);Hand_zBasis_f(i,:)];
     
   %MATRIZ DE ROTACIÓN DE LA ARTICULACIÓN DE LA MUÑECA A PARTIR DEL LAS
   %MATRICES DE ROTACIÓN DE LAS FALANGES DISTAL Y PROXIMAL
   %sin filtro
     R_Wrist(:,:,i)=mtimes(R_Hand(:,:,i),transpose(R_Arm(:,:,i)));
   %Con filtro
     R_Wrist_f(:,:,i)=mtimes(R_Hand_f(:,:,i),transpose(R_Arm_f(:,:,i)));
end

%CALCULO DE LOS ÁNGULOS CARDAN EULER A PARTIR DE LA MATRIZ DE ROTACIÓN MEDIANTE LA FUNCIÓN
%Angulo_Cardan_Euler(R)
%Angulo de la mano Respecto al sistema de referencia global
[Alpha_WG,Beta_WG,Gamma_WG] = Angulo_Cardan_Euler(R_Hand,Ld); 

%Angulo de la muñeca respecto al sistema de referencia local de eslabon proximal y distal
%sin filtro
[Alpha_WL,Beta_WL,Gamma_WL] = Angulo_Cardan_Euler(R_Wrist,Ld);
%con filtro
[Alpha_WL_f,Beta_WL_f,Gamma_WL_f] = Angulo_Cardan_Euler(R_Wrist_f,Ld);

%% Graficas con resultados
figure(6)
plot(Tiempo,rad2deg(-Alpha_WL),'r'),hold on,plot(Tiempo,rad2deg(Beta_WL),'g'),hold off,legend('\alpha','\beta')
xlabel('Tiempo [s]'),ylabel('Angulo [grados]'),title('Angulos Alpha y Beta sin filtro'),grid on;

figure(7)
plot(Tiempo,rad2deg(-Alpha_WL_f),'r'),hold on,plot(Tiempo,rad2deg(Beta_WL_f),'g'),hold off,legend('\alpha','\beta')
xlabel('Tiempo [s]'),ylabel('Angulo [grados]'),title('Angulos {/alpha} y /beta'),grid on;

%% Calculo de la velocidad de los ángulos cardan euler
%Derivadas en  termino de las derivadas de los angulo cardan utilizando la
%eq. 281 del libro Research Methods in biomechanic

%Calculo de la velocidad angular articular de la muñeca
 [vxyz_W] = Velocidad_Angular_CE(-Alpha_WL,Beta_WL,Gamma_WL,Tiempo,Ld);

%Calculo de la velocidad angular articular de la muñeca con filtro
 [vxyz_W_f] = Velocidad_Angular_CE(-Alpha_WL_f,Beta_WL_f,Gamma_WL_f,Tiempo,Ld);  
 
figure(8)
plot(Tiempo,rad2deg(vxyz_W(2,:)),'r'),hold on,plot(Tiempo,rad2deg(vxyz_W_f(2,:)),'c'),hold off,legend('vxyz_W','vxyz_W_f'),
xlabel('Tiempo [s]'),ylabel('Velocidad [grados/s]'),title('Velocidad angular 2D y 3D');grid on;

figure(9)
plot(Tiempo,rad2deg(vxyz_W_f(1,:)),'r'),hold on,plot(Tiempo,rad2deg(vxyz_W_f(2,:)),'g'),hold off,legend('\alpha','\beta'),
xlabel('Tiempo [s]'),ylabel('Velocidad angular [grados/s]'),title('Velocidad angular de la muñeca');grid on;

%% Calculo de la aceleración de la muñeca
%sin filtro
[axyz_W] = Aceleracion_Angular_CE(vxyz_W,Tiempo,Ld);
%con filtro
[axyz_W_f] = Aceleracion_Angular_CE(vxyz_W_f,Tiempo,Ld);

figure(11)
plot(Tiempo,rad2deg(axyz_W(2,:)),'r'),hold on,plot(Tiempo,rad2deg(axyz_W_f(2,:)),'g'),hold off,legend('\alpha_h','\beta_h'),
xlabel('Tiempo [s]'),ylabel('Aceleración angular [grados/s^2]'),title('Aceleración angular de la muñeca'),grid on;

%% dinámica inversa de una junta universal 
%Función que entrega la longitud de la mano, como entrada posición de la
%muñeca (W) y posición final del metacarpio (M_end), falange proximal
%(P_end), falange media (Mi_end) y falange distal (D_end).
[L] = Longitud_Mano(W,M_end,P_end,Mi_end,D_end); %[cm]
P=53; %peso corporal total

%Función que calcula torques aproximados M=[M_alpha, M_beta] y torques exactos M_alpha y M_beta
%si se ingresa true calcula con parámetros dinámicos de hombre y false de mujer
[M_h,M_alpha, M_beta] = Dinamica_Inversa(true,-Alpha_WL_f,Beta_WL_f,vxyz_W_f, axyz_W_f, L, P);
[M_m,M_alpha1,M_beta1] = Dinamica_Inversa(false,-Alpha_WL_f,Beta_WL_f,vxyz_W_f, axyz_W_f, L, P);

%error cuadrático medio normalizado entre torques exactos y aproximados
%NMSE_beta=mean((M_beta-M_h(2,:)).^2)/mean(M_beta.^2)
%NMSE_alpha=mean((M_alpha-M_h(1,:)).^2)/mean(M_alpha.^2)

%% Valores máximos y minimos utilizados para normalizar los datos
min12= -2.3421; max12=1.9251; min34=-1.3443; max34=1.1023; min56=-64.0526; 
max56=55.5066; min7=0.2441; max7=0.4269; min810=-0.2899; max810=2.1599;

%función para crear una base de datos con 10 columnas, las primera 6
%corresponden a las entradas: M_alpha, M_beta, alpha, beta, aceleración
%alpha y beta y las ultimas 4 columnas son las salidas m, K11, K12 y K22.
[Datos_io,~,~] = crear_datos(true,Alpha_WL_f, Beta_WL_f, axyz_W_f, L, P);
%se normalizan los datos
[Datos_n] = normalizar(Datos_io,min12,max12,min34,max34,min56,max56,min7,max7,min810,max810);
%se hace permutación aleatoria las instancias
Datos_n_r=Datos_n(randperm(size(Datos_n, 1)), :);

%% Datos de prueba 
Datos_test=Datos_n_r(700000:800000,:);
Input_test=Datos_test(:,1:6);
%salida para la masa
Output_testm=Datos_test(:,7);
%salida para K11
Output_test11=Datos_test(:,8);
%salida para K12
Output_test12=Datos_test(:,9);
%salida para K22
Output_test22=Datos_test(:,10);


%% Estimación con ANFIS
%Entradas [M_alpha, M_beta, alpha, beta, aceleración alpha, aceleración beta]
%anfisK11=evalfis(fis_K11_hybrid,Input_test(1:50,:));
%anfisK12=evalfis(fis_K12_hybrid,Input_test(1:50,:));
%anfisK22=evalfis(fis_K22_hybrid,Input_test(1:50,:));
%anfism=evalfis(fis_m_hybrid,Input_test(1:50,:));
%% Estimación con redes neuronales
%Entradas [M_alpha, M_beta, alpha, beta, aceleración alpha, aceleración beta]
[Y11,Xf11,Af11] = NN_LM_K11(Input_test(1:50,:));
[Y12,Xf12,Af12] = NN_LM_K12(Input_test(1:50,:));
[Y22,Xf22,Af22] = NN_LM_K22(Input_test(1:50,:));
[Ym,Xm,Am] = NN_LM_m(Input_test(1:50,:));
%% Calculo de R^2 y raíz del error cuadrático medio 
Y0=Y12; %salida de la red
Y=Output_test12(1:50); %salida real

RMSE = sqrt(mean((Y0-Y).^2));
Ybar=mean(Y0);
R=1-(sum((Y0-Y).^2)/sum((Y0-Ybar).^2));

%% Normalización inversa de los datos
%se ingresan las salidas de los modelos y se obtienen sin normalización, lo
%que permite ver los valores en sus unidades originales
[Y11a,Y12a,Y22a,Yma] = Normalizacion_inversa(Y11,Y12,Y22,Ym,min7,max7,min810,max810);
[Output_test11a,Output_test12a,Output_test22a,Output_testma] = Normalizacion_inversa(Output_test11,Output_test12,Output_test22,Output_testm,min7,max7,min810,max810);

subplot(2,2,1)
plot(Y11,'g'),hold on,plot(Output_test11(1:50,:),'c')
legend('K_{\alpha \alpha}^*','K_{\alpha \alpha}')
title('Modelo 1: Parámetro K_{\alpha \alpha}')

subplot(2,2,2)
plot(Y12, 'g'),hold on,plot(Output_test12(1:50,:),'c')
legend('K_{\alpha \beta}^*','K_{\alpha \beta}')
title('Modelo 2: Parámetro K_{\alpha \beta}')

subplot(2,2,3)
plot(Y22,'g'),hold on,plot(Output_test22(1:50,:),'c')
legend('K_{\alpha \beta}^*','K_{\alpha \beta}')
title('Modelo 3: Parámetro K_{\beta \beta}')

subplot(2,2,4)
plot(Ym,'g'),hold on,plot(Output_testm(1:50,:),'c')
legend('m^*','m')
title('Modelo 4: Parámetro m')






