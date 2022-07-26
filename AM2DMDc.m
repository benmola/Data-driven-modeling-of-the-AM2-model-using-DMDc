close all   
clc         
%======================================================================================
m1 = 0.5    ; %   Maximum acidogenic biomass growth rate  
K1 = 7.1    ; %   Half-saturation constant associated with S1
m2 = 0.74   ; %   Maximum methanogenic biomass growth rate
Ki = 80     ; %   Inhibition constant associated with S2
K2 = 9.28   ; %   Half-saturation constant associated with S2
ALPHA = 1   ; %   Proportion of dilution rate for bacteria      
S1in = 20   ; %   S1 Substrate concentration in the input
S2in = 150  ; %   S2 VFA concentration in the input 
k1 = 42.14  ; %   Yield for COD degradation
k2 = 116.5  ; %   Yield for VFA production
k3 = 268    ; %   Yield for VFA consumption
%======================================Simulating the AM2 Model ================================================
KSI_0 = [2;5;1;1]                                                         ;
Tspn = 1:1:121                                                           ;                                  
[Time,KSIP] = ode45(@AM2Model,Tspn,KSI_0)                             ;
AM2_DATA_DMDc = KSIP' ;
% load AM2_DATA_DMDc.mat                               ; 
Ts = 0:1:119                                     ;                                  
Tf = 1:120 ;
% Organizing the data matrices X1 and X2 
X1 = AM2_DATA_DMDc(:, 1:end-1)                      ;
X2 = AM2_DATA_DMDc(:, 2:end)                        ; 
% Control input u = D
D = [0.4*ones(1,40),0.3*ones(1,40),0.2*ones(1,40)]  ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  D  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
figure(1)
lw = 3;
plot(Ts,D,'g','linewidth',lw)                                           ;
legend('Dulution rate','interpreter','latex','linewidth',lw); xlabel('Time [day]','interpreter','latex','linewidth',lw);ylabel('Control input: u = D ($d^{-1}$)' ,'interpreter','latex','linewidth',lw);                                          
axis([0 120 0.19 0.41])                                                   ;
axis square                                                               ;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply the DMDc Algorithm the the set of data X1 and X2 
dt = 1; % one day as a deltha t dt 
[sysmodel_DMDc,U,Up,Ar,Br,G,A_reduced] = DMDc(X1,X2,D,dt)                  ;
xDMDc = lsim(sysmodel_DMDc,D',Ts,X1(:,1))                                ;
X_DMDc = xDMDc'                                                           ; % X_DMDc 
C = [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1]                                     ; % C                    
Y  = C*X_DMDc                                                             ; % Y      
Y1 = Y(1,:)                                                               ; % Y1 
Y2 = Y(2,:)                                                               ; % Y2 
Y3 = Y(3,:)                                                               ; % Y3 
Y4 = Y(4,:)                                                               ; % Y4 
%%

% Plotting the results 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  S1  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)

hold on                                                                   ;
plot(Tf,X1(1,:),'b','LineWidth',lw)                                     ;         
L1 = 'S_1'                                                                ;
plot(Tf,Y1(:,:),':r','linewidth',lw)                                    ;
L2 = 'S_1(DMDc)'                                                          ;                                                                                                  
legend(L1,L2)                                                          ;                                                       
axis square                                                               ;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  X1  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
hold on                                                                   ;
plot(Tf,X1(2,:),'b','LineWidth',lw)                                     ;         
L1 = 'S_1'                                                                ;
plot(Tf,Y2(:,:),':r','linewidth',lw)                                    ;
L2 = 'S_1(DMDc)'                                                          ;                                                                                                  
legend(L1,L2)                                                          ;                                                       
axis square  
%%
                                                            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  S2  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4)
hold on                                                                   ;
plot(Tf,X1(3,:),'b','LineWidth',lw)                                     ;         
L1 = 'S_1'                                                                ;
plot(Tf,Y3(:,:),':r','linewidth',lw)                                    ;
L2 = 'S_1(DMDc)'                                                          ;                                                                                                  
legend(L1,L2)                                                          ;                                                      
axis square  
                                                            

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  X2  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5)
hold on                                                                   ;
plot(Tf,X1(4,:),'b','LineWidth',lw)                                     ;         
L1 = 'S_1'                                                                ;
plot(Tf,Y4(:,:),':r','linewidth',lw)                                    ;
L2 = 'S_1(DMDc)'                                                          ;                                                                                                  
legend(L1,L2)                                                          ;
                                                       
axis square  
                                                             

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  KSI,D  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
figure(6)

                                       

subplot(3,2,3)
hold on                                                                   ;
plot(Tf,X1(1,:),'b','LineWidth',lw)                                     ;         
L1 = '$S_1$' ;                                                                
plot(Tf,Y1(:,:),':r','linewidth',lw)                                    ;
L2 = '$S_1(DMDc)$'; 
legend(L1,L2,'interpreter','latex');
xlabel('Time [day]','interpreter','latex','linewidth',lw);title('Substrate concentration $S_1$ ($g/L$)' ,'interpreter','latex','linewidth',lw);
                                                       

subplot(3,2,4)
hold on                                                                   ;
plot(Tf,X1(2,:),'b','LineWidth',lw)                                     ;         
L1 = '$X_1$' ;                                                               
plot(Tf,Y2(:,:),':r','linewidth',lw)                                    ;
L2 = '$X_1(DMDc)$';                                                                                                                  
legend(L1,L2,'interpreter','latex');
xlabel('Time [day]','interpreter','latex','linewidth',lw);title('Biomass concentration $X_1$ ($g/L$)' ,'interpreter','latex','linewidth',lw);
                                                      

subplot(3,2,5)
hold on                                                                   ;
plot(Tf,X1(3,:),'b','LineWidth',lw)                                     ;         
L1 = '$S_2$' ;                                                               
plot(Tf,Y3(:,:),':r','linewidth',lw)                                    ;
L2 = '$S_2(DMDc)$';                                                                                                                  
legend(L1,L2,'interpreter','latex') ;
xlabel('Time [day]','interpreter','latex','linewidth',lw);title('Substrate concentration $S_2$ ($mmol/L$)' ,'interpreter','latex','linewidth',lw);

                                                      

subplot(3,2,6)
hold on                                                                   ;
plot(Tf,X1(4,:),'b','LineWidth',lw)                                     ;         
L1 = '$X_2$' ;                                                               
plot(Tf,Y4(:,:),':r','linewidth',lw)                                    ;
L2 = '$X_2(DMDc)$' ;                                                                                                                
legend(L1,L2,'interpreter','latex')                             ;
xlabel('Time [day]','interpreter','latex','linewidth',lw);title('Biomass concentration $X_2$ ($g/L$)' ,'interpreter','latex','linewidth',lw);
subplot(3,1,1)
hold on                                                                   ;
plot(Ts,D,'g','linewidth',lw)                                           ;
legend('D : Dillution rate','interpreter','latex') ;
xlabel('Time [day]','interpreter','latex','linewidth',lw);title('Control input: u = D ($d^{-1}$)' ,'interpreter','latex','linewidth',lw);
axis([0 120 0.1 0.41]) 