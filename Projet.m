
%% -------PROJET MODULATION DEMODULATION NUMERIQUE--------------
all clear 
close all
clc

%% ----------Modulation MIA en bande de base---------------
 
%========================== filtre d'émission==============================

K= 100;      %-----------facteur d'echantillonnage
vc =1/(2*K) ;%-----------fréquence de coupure 
n=0:6*K;     %----discrétisation échantillons
N = 500;
n_ny =[-(0:3*K),(1:3*K)]; %----discrétisation échantillons 


h_ny = 2*K*vc*sinc(2*vc*(n_ny)); %---RI RIF non causal
h = 2*K*vc*sinc(2*vc*(n-(6*K)/2));   %---filtre d'émission causal
%---On a multiplié par K, pour amplifier le filtre d'émission, 
%---ainsi avoir un maximum égale à 1.
% Pour appliquer le retard, on a translaté de -3*K

%-------affichage du filtre d'émission---------
% figure() 
%     stem(h) 
%     title("filtre d'émission Causale")
%     xlabel('temps discret')
%     ylabel('amplitude')
%     legend('h[n]')
%     grid("minor")

% ------Affichage du retard de groupe---------
% figure()
%     grpdelay(h,1,N)


%------spectre du signal d'émission---------
[H,~]=freqz(h,1,N); 
% figure()
%     plot(w/(2*pi),log10(abs(H)))
%     title('Spectre du filtre d émission')
%     xlabel('fréquence numérique')
%     ylabel('energie(dB)')
%     grid("minor")


%==========================Séquence d'amplitude============================
long_a = 100;     % -----nombre de symboles 
a = randi([-1,1],1,long_a); %liste aléatoire de taille 1xlong_a de valeur entière comprise entre [-1 1]
% ---On élimine les valeurs nulles en remplaçant les 0 de position i par 1
% si i/2 est pair sinon par -1
 for i = 1:long_a   
     if(a(i)==0 && mod(i,2)==0)  
         a(i)=1;
     elseif (a(i)==0 && mod(i,2)==1)
         a(i)= -1;
     end    
 end

%----------Signal d'Amplitude par interpolation à zéro-----------
Signal_Ampl = zeros(1,long_a*K);
%-------Initialement le signal est nul, on le remplit avec la sequence par
%interpolation 
for i=1:long_a
    Signal_Ampl(i*K) = a(i);
end

%--------------affichage du signal d'amplitude-------------
% figure()
%    stem(Signal_Ampl)
%    title("signal d'amplitude avec 10 symboles ")
%     xlabel('temps discret')
%     ylabel('Amplitude')
%     legend("sequence d'amplitude")
%     grid("minor")



%====================Modulation par filtrage===============================

%--------------Modulation en bande de base----------------
% On sur-échantillonne le signal d'Amplitude pour augmenter la longueur du 
% signal modulé,afin de couvrir tout le signal d'amplitude. ce dernier
%doit être décaler de 300 dû au retard de groupe
% On ajoute 400 échantillons pour combler largement le retard de 800
signal_plus = [Signal_Ampl, zeros(1,800)]; 

retard = 300; 
n_Signal_ampl=retard+1:(length(Signal_Ampl)+retard); % ---sequence d'échantillons 
% pour le signal d'amplitude avec l'application du retard 

signal_module = filter(h,1,signal_plus); % --filtrage du signal d'amplitude sur-échant

%---------Affichage signal de modulation-------------
% figure()
%     plot(signal_module)
%     hold("on")
%     stem(n_Signal_ampl,Signal_Ampl)
%     title("signal de modulation")
%     xlabel('temps discret')
%     ylabel('Amplitude')
%     legend('signal modulé', "signal d'amplitude")
%     grid("minor")


%-----------Spectre du signal de moludation--------------
% [S,~]=freqz(signal_module,1,500);
% figure()
%     plot(w/(2*pi),20*log10(abs(H)))
%     hold("on")
%     plot(w/(2*pi),20*log10(abs(S)))
%     title('Spectre du signal de modulation')
%     xlabel('fréquence numérique')
%     ylabel('energie(dB)')
%     legend("spectre filtre d'émission","spectre signal modulé")
%     grid("minor")


%---------calcul du taux d'erreur-------------------
% On échantillonne le signal modulé par le facteur K avec un retard de 3*K
% + K echantillons: les 3*K corresponds au retard de groupe et le K est dû
% au fait que le signal d'amplitude construit débute toujours avec 100
% valeurs nulles.
% l'échantillonnage est faite que sur la partie qui couvre le signal
% d'amplitude

Echant_s = sign(signal_module(retard+K:K:length(Signal_Ampl)+retard));
error1 = sum(Echant_s ~=a); % calcul d'erreur par comparaison
taux_error_BandeBase = (error1/long_a)*100; % taux d'erreur en pourcentage



% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<FIN PARTIE 1 >>>>>>>>>>>>>>>>>>>
%---------------------------------------------------------------------------------------------------------------------------------------------

%% ------------Modulation MIA sur fréquence porteuse------------
 
vp=0.2;       % --- fréquence porteuse
det_v = 0.05; % ---fréquence de "garde"
n_s = 1:length(signal_module);  % ----discrétisation d'échantillons

s_module_port = signal_module.*cos(2*pi*vp*n_s); % le signal sur fréquence porteuse
[S_port,~]=freqz(s_module_port,1,500); % calcul de gain de s_module_port
 
%------Affichage du signal modulé sur porteuse------------------
% figure()
%     plot(s_module_port)
%     hold("on")
%     stem(n_Signal_ampl,Signal_Ampl)

%------Spectre du signal modulé sur porteuse---------------
% figure()
%     plot(w/(2*pi),20*log10(abs(H)))
%     hold("on")
%     plot(w/(2*pi),20*log10(abs(S_port)))
%     title('Spectre')
%     xlabel('fréquence numérique')
%     ylabel('energie(dB)')
%     legend("impulsion", "signal modulation porteuse")
%     grid("minor")

a_2 = randi([-1,1],1,long_a);  %---construction sequence 2
 for i = 1:long_a
     if(a_2(i)==0 && mod(i,2)==0)
         a_2(i)=-1;
     elseif (a_2(i)==0 && mod(i,2)==1)
         a_2(i)= 1;
     end    
 end

signal_Ampl_2 = zeros(1,6*K);
for i=1:long_a
    signal_Ampl_2(i*K) = a_2(i);
end


a_3 = randi([-1,1],1,long_a); %--- construction séquence 3
 for i = 1:long_a
     if(a_3(i)==0 && mod(i,2)==0)
         a_3(i)=1;
     elseif (a_3(i)==0 && mod(i,2)==1)
         a_3(i)= -1;
     end    
 end

signal_Ampl_3 = zeros(1,6*K);
for i=1:long_a
    signal_Ampl_3(i*K) = a_3(i);
end

%sur-échantillonnage des signaux d'amplitude
seq_2 = [signal_Ampl_2, zeros(1,800)];
seq_3 = [signal_Ampl_3, zeros(1,800)];

% codage des signaux en bande de base
s_2 = filter(h,1,seq_2); 
s_3 = filter(h,1,seq_3);

%modulation sur fréquence porteuse avec garde
signal_garde_1 = s_2.*cos(2*pi*(vp+det_v)*n_s);
signal_garde_2 = s_3.*cos(2*pi*(vp-det_v)*n_s);

% somme des 3 signaux porteurs
signal_somme = s_module_port +signal_garde_1+ signal_garde_2;

% ----------calcul de gain-------------
[S_G_1,~]=freqz(signal_garde_2,1,500);
[S_G_2,~]=freqz(signal_garde_1,1,500);
[S_som,~]=freqz(signal_somme,1,500);

%--------Affichage des spectres séparement -------------
% figure()
%     hold on
%     plot(w/(2*pi),mag2db(abs(H)));
%     plot(w/(2*pi),mag2db(abs(S_port)));
%     plot(w/(2*pi),mag2db(abs(S_G_2)));
%     plot(w/(2*pi),mag2db(abs(S_G_1)));
%     title('Spectre')
%     xlabel('fréquence numérique')
%     ylabel('energie(dB)')
%     legend("impulsion", "porteuse","νp + δν", "νp - δν")
%     grid("minor")
%     hold off


% ---------Affichage du spectre porteuse avec les "gardes"-------
% figure()
%     hold on
%     plot(w/(2*pi),mag2db(abs(H)));
%     plot(w/(2*pi),mag2db(abs(S_som)));
%     title('Spectre')
%     xlabel('fréquence numérique')
%     ylabel('energie(dB)')
%     legend("impulsion", "porteuse + 'gardes'")
%     grid("minor")
%     hold off

 
%=========================Démodulation===================================

% --------------Rémultiplier le signal central par la porteur-----------

s_a_demod = (1/K)*signal_module.*(cos(2*pi*vp*n_s)).^2; 

[S_a_demod,~]=freqz(s_a_demod,1,500);

%-----spectre du signal à démoduler-----------
% figure()
%     plot(w/(2*pi),20*log10(abs(H)))
%     hold("on")
%     plot(w/(2*pi),20*log10(abs(S_a_demod)))
%     title('Spectre')
%     xlabel('fréquence numérique')
%     ylabel('energie(dB)')
%     legend("impultion", "démodulation cohérente")
%     grid("minor")


%-------------- Filtre pass-bas de demodulation--------------
% si on peut bien isoler le signal porteur, on prend un filtre de passe
% bande la lagueur de la bande de garde

v0 = 2*det_v;  %vc ------fréquence de coupure du filtre passe-bas
ordre = 2;  % ordre du filtre

% filtre passe_bas avec fênetre de Hann() pour diminuer l'amplitude des
% lobes et éviter une variation de l'amplitude du signal démodulé 
h_passBas = fir1(ordre,v0,"low",hann(ordre+1)); 


[H_passBas,~]=freqz(h_passBas,1,500);

%------------spectre du filtre-------------------
% figure()
%     plot(w/(2*pi),20*log10(abs(H)))
%     hold("on")
%     plot(w/(2*pi),20*log10(abs(S_a_demod)))
%     plot(w/(2*pi),20*log10(abs(H_passBas)),'g')
%     title('Spectre')
%     xlabel('fréquence numérique')
%     ylabel('energie(dB)')
%     legend("implusion", "demodulation","passeBas+hann()")
%     grid('minor')

%---------- filtrage pour démodulation----------------
% signal démodulé, on amplifie le signal pour observer et vérifier la condition de
% Nyquist (le filtre n'a pas de partie amplificateur)
signal_demodule = 2*K*filter(h_passBas,1,s_a_demod);

[g_2,~]=grpdelay(h_passBas,1,500); %calcul de retard de groupe
retard_2 = mean(g_2) + retard;
%Compensation du retard de groupe pour le signal en bande de base
n_2=retard_2+1:(length(Signal_Ampl)+retard_2);
[S_demo,~]=freqz(signal_demodule,1,500);

% ---------Affichage du signal de démodulation-----------------
% figure()
%     plot(signal_demodule)    
%     hold("on")
%     stem(n_2,Signal_Ampl)
%     title('signal de démodulation')
%     xlabel('fréquence numérique')
%     ylabel('energie(dB)')
%     legend(" signal démodulé", "signal d'Amplitude")
%     grid('minor')
    

%-----------Calcul taux d'erreur---------------------------
Echant_s_d = sign(signal_demodule(retard_2+K:K:length(Signal_Ampl)+retard_2));
error2_demod = sum(Echant_s_d ~=a);
%error2_demod = sum(round(signal_module) ~= round(signal_demodule));
taux_error_demod = (error2_demod/long_a)*100;

%avec une ordre=2, la démodulation ne presente pas d'erreur, si on compense
%le décalage dû au retard de groupe. 
%Mais le signal de bande de bande ne se superpose pas sur le signal de
%démodulation
% + l'ordre est élevé + le retard est important dû au temps de calcul


%<<<<<<<<<<<<<<<<<<<<<<<<<FIN PARTIE 2 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%--------------------------------------------------------------------------------------------------------------------------------------------

%% ------------------Etude de la resistance au bruit----------------------
 
%------------------signal bruité--------
 o =  5; %Amplitude du bruit
 %--------- signal de démodulation bruité avec un bruit gaussien--------
 signal_bruite = s_a_demod + o*randn(1,length(s_a_demod)); 

%-----------filtre de démodulation -----------
v02 = 2*det_v;
ordre2 = 10;
h_passBas2 = fir1(ordre2,v02,"low",hann(ordre2+1));
[H_passBas2,w]=freqz(h_passBas2,1,500);

%---------retard de groupe-------------
[g_3,w_3]=grpdelay(h_passBas2,1,500);% calcul de retard de groupe
retard_3 = mean(g_3); % le retard de groupe
n_3=retard_3+1:(length(signal_module)+ retard_3);%temps discrét pour le signal de bande de base
n_3_2=retard_3+1:(length(Signal_Ampl)+ retard_3); %temps discrét pour le signal d'amplitude

%----------filtrage pour démodumation------------------
signal_bruite_demo = 2*K*filter(h_passBas2,1,signal_bruite);
[Signal_bruite_demo,~]=freqz(signal_bruite_demo,1,500);



%------------Spectre du filtre et du signal bruité-------------
% figure()
%     plot(w/(2*pi),20*log10(abs(H)))
%     hold("on")
%     plot(w/(2*pi),20*log10(abs(Signal_bruite_demo)))
%     plot(w/(2*pi),20*log10(abs(H_passBas2)),'g')
%     title('Spectre')
%     xlabel('fréquence numérique')
%     ylabel('energie(dB)')
%     legend("impulsion", "demodulation bruité", "passeBas")
%     grid('minor')


%---------Affichage du signal bruité-----------    
% figure()
%     plot(signal_bruite_demo)    
%     hold("on")
%     plot(n_3,signal_module)
%     stem(n_3_2,Signal_Ampl,'b') 
%     title('signal de démodulation bruité')
%     xlabel('fréquence numérique')
%     ylabel('energie(dB)')
%     legend(" signal bruité", "signal bande de base","signal d'amplitude")
%     grid('minor')


  %------------- calcul erreur -------------
%taux_error_demod_bruit = erreur_demod_bruit/length(signal_module);
Echant_demod_bruit= sign(signal_bruite_demo(retard_3+K:K:length(Signal_Ampl)+retard_3));
error3_demod_bruit = sum(Echant_s_d ~=a);
taux_error_demod_bruit = (error3_demod_bruit/long_a)*100;

    

