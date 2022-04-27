function plot_hand(Input,Frame_start,Step,Frame_end)
%Función que permite crear los vectores que definen las falanges dista,
%media, proximal y metacarpianos, tomar la posición de la muñeca y del
%codo respecto al sistema de coordenadas cartesiano del sensor leap motion

Datos = table2array(Input);
    %% Se definen matrices de (2,3,numero de datos), correspondientes a las Falanges que se desean graficar

%Falanges del dedo meñique

Ld=length(Datos(:,1));
P_fd=zeros(2,3,Ld);
P_fm=zeros(2,3,Ld);
P_fp=zeros(2,3,Ld);
P_mt=zeros(2,3,Ld);

%Falanges del dedo anular
R_fd=zeros(2,3,Ld);
R_fm=zeros(2,3,Ld);
R_fp=zeros(2,3,Ld);
R_mt=zeros(2,3,Ld);

%Falanges del dedo medio
M_fd=zeros(2,3,Ld);
M_fm=zeros(2,3,Ld);
M_fp=zeros(2,3,Ld);
M_mt=zeros(2,3,Ld);

%Falanges del dedo índice
I_fd=zeros(2,3,Ld);
I_fm=zeros(2,3,Ld);
I_fp=zeros(2,3,Ld);
I_mt=zeros(2,3,Ld);

%Falanges del dedo pulgar
T_fd=zeros(2,3,Ld);
T_fm=zeros(2,3,Ld);
T_fp=zeros(2,3,Ld);
T_mt=zeros(2,3,Ld);

% Falange que describe el antebrazo
forearm=zeros(2,3,Ld);
hand=zeros(2,3,Ld);

%líneas que unen los nudillos de los dedos índice, medio, anular y meñique
A=zeros(2,3,Ld);
B=zeros(2,3,Ld);
C=zeros(2,3,Ld);

%líneas que unen el inicio de los metacarpiano de los dedos índice, medio,
%anular y meñique, en la zona de la muñeca
Aw=zeros(2,3,Ld);
Bw=zeros(2,3,Ld);
Cw=zeros(2,3,Ld);

Hand_Basisx=zeros(2,3,Ld);
Hand_Basisy=zeros(2,3,Ld);
Hand_Basisz=zeros(2,3,Ld);

Arm_Basisx=zeros(2,3,Ld);
Arm_Basisy=zeros(2,3,Ld);
Arm_Basisz=zeros(2,3,Ld);

M_Basisx=zeros(2,3,Ld);
M_Basisy=zeros(2,3,Ld);
M_Basisz=zeros(2,3,Ld);

LM_Basisx=zeros(2,3,Ld);
LM_Basisy=zeros(2,3,Ld);
LM_Basisz=zeros(2,3,Ld);


Arm_0=zeros(Ld,3);

W=table2array(Input(:,{'wrist_position_x','wrist_position_y','wrist_position_z'}));
  % Ubicación de las posición (en x, y, z ) del codo
E=table2array(Input(:,{'elbow_position_x','elbow_position_y','elbow_position_z'}));
for k=1:Ld
    Arm_0(k,:)=E(k,:)+(W(k,:)-E(k,:))*0.5;
end




%% Creación los puntos de las articulaciones a graficar
for i=Frame_start:Step:Frame_end
    %% Dedo meñique
    % Ubicación de las posición final (en x, y, z ) de las articulaciones de las falanges 
    %distal (afd),media (afm),  proximal (afp), metacarpianos (amt)
    %P_afd=Datos(i,77:79);
    P_afd=table2array(Input(i,{'Pinky_fin_dist_end_x','Pinky_fin_dist_end_y','Pinky_fin_dist_end_z'}));
    P_afm=table2array(Input(i,{'Pinky_fin_inter_end_x','Pinky_fin_inter_end_y','Pinky_fin_inter_end_z'}));
    P_afp=table2array(Input(i,{'Pinky_fin_prox_end_x','Pinky_fin_prox_end_y','Pinky_fin_prox_end_z'}));
    P_amt=table2array(Input(i,{'Pinky_fin_meta_end_x','Pinky_fin_meta_end_y','Pinky_fin_meta_end_z'}));
    P_amtw=table2array(Input(i,{'Pinky_fin_meta_start_x','Pinky_fin_meta_start_y','Pinky_fin_meta_start_z'}));

    % Se crean las falanges uniendo los puntos de las articulaciones
    % correspondientes 
    P_fd(:,:,i)=[P_afd;P_afm];
    P_fm(:,:,i)=[P_afm;P_afp];
    P_fp(:,:,i)=[P_afp;P_amt];
    P_mt(:,:,i)=[P_amt;P_amtw];

    %% Dedo anular
    % Ubicación de las posición final (en x, y, z ) de las articulaciones de las falanges 
    %distal (afd),media (afm),  proximal (afp), metacarpianos (amt)
    R_afd=table2array(Input(i,{'Ring_fin_dist_end_x','Ring_fin_dist_end_y','Ring_fin_dist_end_z'}));
    R_afm=table2array(Input(i,{'Ring_fin_inter_end_x','Ring_fin_inter_end_y','Ring_fin_inter_end_z'}));
    R_afp=table2array(Input(i,{'Ring_fin_prox_end_x','Ring_fin_prox_end_y','Ring_fin_prox_end_z'}));
    R_amt=table2array(Input(i,{'Ring_fin_meta_end_x','Ring_fin_meta_end_y','Ring_fin_meta_end_z'}));
    R_amtw=table2array(Input(i,{'Ring_fin_meta_start_x','Ring_fin_meta_start_y','Ring_fin_meta_start_z'}));


    % Se crean las falanges uniendo los puntos de las articulaciones
    % correspondientes 
    R_fd(:,:,i)=[R_afd;R_afm];
    R_fm(:,:,i)=[R_afm;R_afp];
    R_fp(:,:,i)=[R_afp;R_amt];
    R_mt(:,:,i)=[R_amt;R_amtw];
    
    %% Dedo medio
    %Ubicación de las posición final (en x, y, z ) de las articulaciones de las falanges 
    %distal (afd),media (afm),  proximal (afp), metacarpianos (amt)
     
    M_afd=table2array(Input(i,{'Middle_fin_dist_end_x','Middle_fin_dist_end_y','Middle_fin_dist_end_z'}));
    M_afm=table2array(Input(i,{'Middle_fin_inter_end_x','Middle_fin_inter_end_y','Middle_fin_inter_end_z'}));
    M_afp=table2array(Input(i,{'Middle_fin_prox_end_x','Middle_fin_prox_end_y','Middle_fin_prox_end_z'}));
    M_amt=table2array(Input(i,{'Middle_fin_meta_end_x','Middle_fin_meta_end_y','Middle_fin_meta_end_z'}));
    M_amtw=table2array(Input(i,{'Middle_fin_meta_start_x','Middle_fin_meta_start_y','Middle_fin_meta_start_z'}));
    
    %Se crean las falanges uniendo los puntos de las articulaciones correspondientes 
    M_fd(:,:,i)=[M_afd;M_afm];
    M_fm(:,:,i)=[M_afm;M_afp];
    M_fp(:,:,i)=[M_afp;M_amt];
    M_mt(:,:,i)=[M_amt;M_amtw];

    %% Dedo índice
    % Ubicación de las posición final (en x, y, z ) de las articulaciones de las falanges 
    %distal (afd),media (afm),  proximal (afp), metacarpianos (amt)
    
    I_afd=table2array(Input(i,{'Index_fin_dist_end_x','Index_fin_dist_end_y','Index_fin_dist_end_z'}));
    I_afm=table2array(Input(i,{'Index_fin_inter_end_x','Index_fin_inter_end_y','Index_fin_inter_end_z'}));
    I_afp=table2array(Input(i,{'Index_fin_prox_end_x','Index_fin_prox_end_y','Index_fin_prox_end_z'}));
    I_amt=table2array(Input(i,{'Index_fin_meta_end_x','Index_fin_meta_end_y','Index_fin_meta_end_z'}));
    I_amtw=table2array(Input(i,{'Index_fin_meta_start_x','Index_fin_meta_start_y','Index_fin_meta_start_z'}));

    % Se crean las falanges uniendo los puntos de las articulaciones
    % correspondientes 
    I_fd(:,:,i)=[I_afd;I_afm];
    I_fm(:,:,i)=[I_afm;I_afp];
    I_fp(:,:,i)=[I_afp;I_amt];
    I_mt(:,:,i)=[I_amt;I_amtw];

    %% Dedo pulgar
    % Ubicación de las posición final (en x, y, z ) de las articulaciones de las falanges 
    %distal (afd),media (afm),  proximal (afp), metacarpianos (amt)
    
    T_afd=table2array(Input(i,{'Thumb_fin_dist_end_x','Thumb_fin_dist_end_y','Thumb_fin_dist_end_z'}));
    T_afm=table2array(Input(i,{'Thumb_fin_inter_end_x','Thumb_fin_inter_end_y','Thumb_fin_inter_end_z'}));
    T_afp=table2array(Input(i,{'Thumb_fin_prox_end_x','Thumb_fin_prox_end_y','Thumb_fin_prox_end_z'}));
    T_amt=table2array(Input(i,{'Thumb_fin_meta_end_x','Thumb_fin_meta_end_y','Thumb_fin_meta_end_z'}));
    T_amtw=table2array(Input(i,{'Thumb_fin_meta_start_x','Thumb_fin_meta_start_y','Thumb_fin_meta_start_z'}));

    % Se crean las falanges uniendo los puntos de las articulaciones
    % correspondientes 
    T_fd(:,:,i)=[T_afd;T_afm];
    T_fm(:,:,i)=[T_afm;T_afp];
    T_fp(:,:,i)=[T_afp;T_amt];
    T_mt(:,:,i)=[T_amt;T_amtw];

    %% Antebrazo que va desde el codo a la muñeca
    % Ubicación de las posición (en x, y, z ) de la muñeca
    W=table2array(Input(i,{'wrist_position_x','wrist_position_y','wrist_position_z'}));
      % Ubicación de las posición (en x, y, z ) del codo
    E=table2array(Input(i,{'elbow_position_x','elbow_position_y','elbow_position_z'})); 
     
   
    % Se crea la falange que describre el antebrazo
    forearm(:,:,i)=[W;E];
    hand(:,:,i)=[M_afd;W];
    

    %% conformación de la mano
    % líneas que unen los nudillos de los dedos
    A(:,:,i)=[P_amt;R_amt];
    B(:,:,i)=[R_amt;M_amt];
    C(:,:,i)=[M_amt;I_amt];
    % líneas que unen el inicio de los metacaparpianos, en la zona de la
    % muñeca
    Aw(:,:,i)=[P_amtw;R_amtw];
    Bw(:,:,i)=[R_amtw;M_amtw];
    Cw(:,:,i)=[M_amtw;I_amtw];

    %% Vectores base
    Hand_0=table2array(Input(i,{'hand_position_x','hand_position_y','hand_position_z'})); 
    
    M_xBasis=table2array(Input(i,{'Middle_x_basis_x','Middle_x_basis_y','Middle_x_basis_z'}));
    M_yBasis=table2array(Input(i,{'Middle_y_basis_x','Middle_y_basis_y','Middle_y_basis_z'}));
    M_zBasis=table2array(Input(i,{'Middle_z_basis_x','Middle_z_basis_y','Middle_z_basis_z'}));

    Arm_xBasis=table2array(Input(i,{'arm_basis_x_x','arm_basis_x_y','arm_basis_x_z'}));
    Arm_yBasis=table2array(Input(i,{'arm_basis_y_x','arm_basis_y_y','arm_basis_y_z'}));
    Arm_zBasis=table2array(Input(i,{'arm_basis_z_x','arm_basis_z_y','arm_basis_z_z'}));

    Hand_xBasis=table2array(Input(i,{'Hand_basis_x_x','Hand_basis_x_y','Hand_basis_x_z'}));
    Hand_yBasis=table2array(Input(i,{'Hand_basis_y_x','Hand_basis_y_y','Hand_basis_y_z'}));
    Hand_zBasis=table2array(Input(i,{'Hand_basis_z_x','Hand_basis_z_y','Hand_basis_z_z'}));
    
    Origin=table2array(Input(i,{'Middle_origin_x','Middle_origin_y','Middle_origin_z'}));
    
   
    
    M_Basisx(:,:,i)=[M_afp;M_afp+M_xBasis*50];
    M_Basisy(:,:,i)=[M_afp;M_afp+M_yBasis*50];
    M_Basisz(:,:,i)=[M_afp;M_afp+M_zBasis*50];
    
    Hand_Basisx(:,:,i)=[Hand_0;Hand_0+Hand_xBasis*50];
    Hand_Basisy(:,:,i)=[Hand_0;Hand_0+Hand_yBasis*50];
    Hand_Basisz(:,:,i)=[Hand_0;Hand_0+Hand_zBasis*50];
    
    Arm_Basisx(:,:,i)=[Arm_0(i,:);(Arm_0(i,:)+Arm_xBasis*50)];
    Arm_Basisy(:,:,i)=[Arm_0(i,:);(Arm_0(i,:)+Arm_yBasis*50)];
    Arm_Basisz(:,:,i)=[Arm_0(i,:);(Arm_0(i,:)+Arm_zBasis*50)];
    
    LM_Basisx(:,:,i)=[[0 0 0];[50 0 0]];
    LM_Basisy(:,:,i)=[[0 0 0];[0 50 0]];
    LM_Basisz(:,:,i)=[[0 0 0];[0 0 50]];
    
    


    %% Grafica
        
    plot3(P_fd(:,1,i),P_fd(:,2,i),P_fd(:,3,i),'-go','LineWidth',3,'MarkerSize',4,'MarkerEdgeColor','b','MarkerFaceColor','b')
    hold on
    plot3(P_fm(:,1,i),P_fm(:,2,i),P_fm(:,3,i),'-go','LineWidth',3,'MarkerSize',4,'MarkerEdgeColor','b','MarkerFaceColor','b')
    hold on
    plot3(P_fp(:,1,i),P_fp(:,2,i),P_fp(:,3,i),'-go','LineWidth',3,'MarkerSize',4,'MarkerEdgeColor','b','MarkerFaceColor','b')
    hold on
    plot3(P_mt(:,1,i),P_mt(:,2,i),P_mt(:,3,i),'-go','LineWidth',3,'MarkerSize',4,'MarkerEdgeColor','b','MarkerFaceColor','b')
    hold on

    plot3(R_fd(:,1,i),R_fd(:,2,i),R_fd(:,3,i),'-go','LineWidth',3,'MarkerSize',4,'MarkerEdgeColor','b','MarkerFaceColor','b')
    hold on
    plot3(R_fm(:,1,i),R_fm(:,2,i),R_fm(:,3,i),'-go','LineWidth',3,'MarkerSize',4,'MarkerEdgeColor','b','MarkerFaceColor','b')
    hold on
    plot3(R_fp(:,1,i),R_fp(:,2,i),R_fp(:,3,i),'-go','LineWidth',3,'MarkerSize',4,'MarkerEdgeColor','b','MarkerFaceColor','b')
    hold on
    plot3(R_mt(:,1,i),R_mt(:,2,i),R_mt(:,3,i),'-go','LineWidth',3,'MarkerSize',4,'MarkerEdgeColor','b','MarkerFaceColor','b')
    hold on

    plot3(M_fd(:,1,i),M_fd(:,2,i),M_fd(:,3,i),'-go','LineWidth',3,'MarkerSize',4,'MarkerEdgeColor','b','MarkerFaceColor','b')
    hold on
    plot3(M_fm(:,1,i),M_fm(:,2,i),M_fm(:,3,i),'-go','LineWidth',3,'MarkerSize',4,'MarkerEdgeColor','b','MarkerFaceColor','b')
    hold on
    plot3(M_fp(:,1,i),M_fp(:,2,i),M_fp(:,3,i),'-go','LineWidth',3,'MarkerSize',4,'MarkerEdgeColor','b','MarkerFaceColor','b')
    hold on
    plot3(M_mt(:,1,i),M_mt(:,2,i),M_mt(:,3,i),'-go','LineWidth',3,'MarkerSize',4,'MarkerEdgeColor','b','MarkerFaceColor','b')
    hold on

    plot3(I_fd(:,1,i),I_fd(:,2,i),I_fd(:,3,i),'-go','LineWidth',3,'MarkerSize',4,'MarkerEdgeColor','b','MarkerFaceColor','b')
    hold on
    plot3(I_fm(:,1,i),I_fm(:,2,i),I_fm(:,3,i),'-go','LineWidth',3,'MarkerSize',4,'MarkerEdgeColor','b','MarkerFaceColor','b')
    hold on
    plot3(I_fp(:,1,i),I_fp(:,2,i),I_fp(:,3,i),'-go','LineWidth',3,'MarkerSize',4,'MarkerEdgeColor','b','MarkerFaceColor','b')
    hold on
    plot3(I_mt(:,1,i),I_mt(:,2,i),I_mt(:,3,i),'-go','LineWidth',3,'MarkerSize',4,'MarkerEdgeColor','b','MarkerFaceColor','b')
    hold on

    plot3(T_fd(:,1,i),T_fd(:,2,i),T_fd(:,3,i),'-go','LineWidth',3,'MarkerSize',4,'MarkerEdgeColor','b','MarkerFaceColor','b')
    hold on
    plot3(T_fm(:,1,i),T_fm(:,2,i),T_fm(:,3,i),'-go','LineWidth',3,'MarkerSize',4,'MarkerEdgeColor','b','MarkerFaceColor','b')
    hold on
    plot3(T_fp(:,1,i),T_fp(:,2,i),T_fp(:,3,i),'-go','LineWidth',3,'MarkerSize',4,'MarkerEdgeColor','b','MarkerFaceColor','b')
    hold on
    plot3(T_mt(:,1,i),T_mt(:,2,i),T_mt(:,3,i),'-go','LineWidth',3,'MarkerSize',4,'MarkerEdgeColor','b','MarkerFaceColor','b')
    hold on

   
    
    plot3(C(:,1,i),C(:,2,i),C(:,3,i),'-go','LineWidth',3,'MarkerSize',4,'MarkerEdgeColor','b','MarkerFaceColor','b')
    hold on
    plot3(Aw(:,1,i),Aw(:,2,i),Aw(:,3,i),'-go','LineWidth',3,'MarkerSize',4,'MarkerEdgeColor','b','MarkerFaceColor','b')
    hold on
    plot3(Bw(:,1,i),Bw(:,2,i),Bw(:,3,i),'-go','LineWidth',3,'MarkerSize',4,'MarkerEdgeColor','b','MarkerFaceColor','b')
    hold on
    plot3(Cw(:,1,i),Cw(:,2,i),Cw(:,3,i),'-go','LineWidth',3,'MarkerSize',4,'MarkerEdgeColor','b','MarkerFaceColor','b')
    hold on
    plot3(forearm(:,1,i),forearm(:,2,i),forearm(:,3,i),'-go','LineWidth',3,'MarkerSize',4,'MarkerEdgeColor','b','MarkerFaceColor','b')
    %hold on
    %plot3(hand(:,1,i),hand(:,2,i),hand(:,3,i),'-ro','LineWidth',3,'MarkerSize',2,'MarkerEdgeColor','r','MarkerFaceColor','r')
    
    
    hold on
    plot3(A(:,1,i),A(:,2,i),A(:,3,i),'-go','LineWidth',3,'MarkerSize',4,'MarkerEdgeColor','b','MarkerFaceColor','b')
    hold on
    plot3(B(:,1,i),B(:,2,i),B(:,3,i),'-go','LineWidth',3,'MarkerSize',4,'MarkerEdgeColor','b','MarkerFaceColor','b')
    hold on
    
    
    
    
    plot3(Hand_Basisx(:,1,i),Hand_Basisx(:,2,i),Hand_Basisx(:,3,i),'-y>','LineWidth',2,'MarkerSize',2,'MarkerEdgeColor','y','MarkerFaceColor','y')
    hold on
    plot3(Hand_Basisy(:,1,i),Hand_Basisy(:,2,i),Hand_Basisy(:,3,i),'-b>','LineWidth',2,'MarkerSize',2,'MarkerEdgeColor','b','MarkerFaceColor','b')
    hold on
    plot3(Hand_Basisz(:,1,i),Hand_Basisz(:,2,i),Hand_Basisz(:,3,i),'-r>','LineWidth',2,'MarkerSize',2,'MarkerEdgeColor','r','MarkerFaceColor','r')
    hold on
    plot3(Arm_Basisx(:,1,i),Arm_Basisx(:,2,i),Arm_Basisx(:,3,i),'-y>','LineWidth',2,'MarkerSize',2,'MarkerEdgeColor','y','MarkerFaceColor','y')
    hold on
    plot3(Arm_Basisy(:,1,i),Arm_Basisy(:,2,i),Arm_Basisy(:,3,i),'-b>','LineWidth',2,'MarkerSize',2,'MarkerEdgeColor','b','MarkerFaceColor','b')
    hold on
    plot3(Arm_Basisz(:,1,i),Arm_Basisz(:,2,i),Arm_Basisz(:,3,i),'-r>','LineWidth',2,'MarkerSize',2,'MarkerEdgeColor','r','MarkerFaceColor','r')
    
    plot3(LM_Basisx(:,1,i),LM_Basisx(:,2,i),LM_Basisx(:,3,i),'-y>','LineWidth',2,'MarkerSize',2,'MarkerEdgeColor','y','MarkerFaceColor','y')
    hold on
    plot3(LM_Basisy(:,1,i),LM_Basisy(:,2,i),LM_Basisy(:,3,i),'-b>','LineWidth',2,'MarkerSize',2,'MarkerEdgeColor','b','MarkerFaceColor','b')
    hold on
    plot3(LM_Basisz(:,1,i),LM_Basisz(:,2,i),LM_Basisz(:,3,i),'-r>','LineWidth',2,'MarkerSize',2,'MarkerEdgeColor','r','MarkerFaceColor','r')
    
    hold off 
    grid on
    title('Gráfica del modelo simplificado de la mano')
    axis equal
    xlabel('X [mm]')
    ylabel('Y [mm]')
    zlabel('Z [mm]')
    drawnow

end



end

