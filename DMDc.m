function [sysmodel_DMDc,U,Up,Ar,Br,G,A_reduced,sys_reduced ] = DMDc(X1,X2,U,dt)

numOutputs = size(X1,1)                                          ;
numInputs = size(U,1)                                            ;
[l,c] = size(U)                                                  ;
red = size(X1,1) - 1                                             ;
Omega = [X1; U]                                                  ;
[U,S,V] = svd(Omega,'econ')                                      ; 
G = X2*V*S^(-1)*U';


[Up,Sp,Vp] = svd(X2,'econ')                                      ; 


Ar = G(:,1:numOutputs)                                           ;
Br = G(:,numOutputs+1:end)                                       ;
Cr = eye(numOutputs)                                             ; 


sysmodel_DMDc = ss(Ar,Br,Cr,zeros(numOutputs,numInputs),1)       ;

G_reduced = Up'*X2*V*S^(-1)*U'                                   ;
A_reduced = G_reduced(:,1:numOutputs)*Up                         ;
A_reduced = A_reduced(1:red,1:red)                               ;

B_reduced = Up'*G(:,numOutputs+1:end)                            ;
B_reduced = B_reduced(1:red,1)                                   ;
C_reduced = eye(red)                                             ;
D_reduced = zeros(3,1)                                           ;
sys_reduced = ss(A_reduced,B_reduced,C_reduced,D_reduced,dt)     ;
