function KSIP = AM2Model(TEMPS,KSI)

% Varriable control input u = D
if TEMPS     <= 40
    D = 0.4         ;
elseif TEMPS <= 80
    D = 0.3         ;
elseif TEMPS <= 120
    D = 0.2         ;
else
    D = 0.2         ;
end

%============================AM2 Model parameters ==========================================================

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

%======================================================================================


%%%%%%%%%%%%%%%%%%%% Monod growth kinetics %%%%%%%%%%%%%%%%%%%%%%%%%%
MU1 = m1 .* ( KSI(1)./(KSI(1) + K1) )                            ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%  Haldane growth kinetics %%%%%%%%%%%%%%%%%%%%%%%%%
MU2 = m2 .* (KSI(3)./((KSI(3).^2/Ki) +KSI(3)+ K2) )              ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

KSIP = zeros(4,1)                                                ;              
KSIP(1) = (D*(S1in - KSI(1))) - (k1*MU1*KSI(2))                  ;
KSIP(2) = (MU1 - ALPHA*D)*KSI(2)                                 ;
KSIP(3) = (D*(S2in-KSI(3)) )+ (k2*MU1*KSI(2)) - (k3*MU2*KSI(4))  ;
KSIP(4) = (MU2 - D)*KSI(4)                                       ;  

end