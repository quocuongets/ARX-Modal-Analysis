%Gerer les donnees temps
X=Fy_mod;   % grinding force - input signal
Y=Ay_mod;   % acceleration - output signal

[Ly,ny]=size(Y);
[Lu,nu]=size(X);
Ts=1/512;
d=ny;

%Structure de donnee
Z=iddata(Y,X,Ts);

for p=2:1:100
%Ordre du modele
% na=10;
% nb=10;
% nk=1;

na=p*ones(ny,ny);
nb=p*ones(ny,nu);
nk=1*ones(ny,nu);

%Modele arx
M=arx(Z,[na nb nk]);

% Partie autoregressive
A=M.A; %Une matrice des parametres

%Matrice des parametres autoregressives
S=[];
for i=1:p
    S=[S -A(:,i+1)];
end

%Matrice d'etat
A1=[S; eye((p-1)*d) zeros((p-1)*d,d)];

% A1=[S; eye((p-d)*d) zeros((p-d),d)];

%Decomposition valeurs propres
[vect,val]=eig(A1);

%Calculer les frequences et taux d'amorts
u=diag(val);
lamda=log(u)/Ts;
omega=zeros(1,length(lamda));
freq=zeros(1,length(lamda));
ksil=zeros(1,length(lamda));
for i=1:length(lamda)
    omega(i)=sqrt(real(lamda(i))^2+imag(lamda(i))^2);
    freq(i)=omega(i)/2/pi;
    ksil(i)=100*abs(real(lamda(i)))/omega(i);
end

%Calculer les modes
amplitude=abs(vect);
for i=1:size(amplitude,2)
    amplitude(:,i)=amplitude(:,i)/amplitude(1,i);
end
mode=amplitude(1:d,:);

figure(1)
plot(freq, p*ones,'.','LineWidth',1.5)
% title('Frequency stabilization diagrams.')
ylabel('Model order')
xlabel('Frequency (Hz)')
hold on

% figure(2)
% for j=1:length(freq)
%     if 0 <= freq(j) <= 100
%     plot(ksil(j),p,'k+')
%     title('Damping ratio stabilization diagrams')
% ylabel('Model order')
% xlabel('Damping ratio')
%     end
% end
% hold on


end

return

% %At order max
% %Retrieve the excitation parameters
% B=M.B;

% %Matrice des parametres eXegenous
% T=[];
% for i=1:p
%     X=[T B(:,:,i+1)];
% end

%Freqresp function
f=[0:1/4:120];
W=2*pi*f;
H=freqresp(M,W);

% Hxy: x output, y input

Hxx=[];
for i=1:size(H,3)
    Hxx(i)=H(1,1,i);
end

Hxy=[];
for i=1:size(H,3)
    Hxy(i)=H(1,2,i);
end

Hxz=[];
for i=1:size(H,3)
    Hxz(i)=H(1,3,i);
end


Hyx=[];
for i=1:size(H,3)
    Hyx(i)=H(2,1,i);
end

Hyy=[];
for i=1:size(H,3)
    Hyy(i)=H(2,2,i);
end

Hyz=[];
for i=1:size(H,3)
    Hyz(i)=H(2,3,i);
end

Hzx=[];
for i=1:size(H,3)
    Hzx(i)=H(3,1,i);
end

Hzy=[];
for i=1:size(H,3)
    Hzy(i)=H(3,2,i);
end

Hzz=[];
for i=1:size(H,3)
    Hzz(i)=H(3,3,i);
end

figure(3)
subplot(2,1,1)
plot(f,abs(Hxx))
ylabel('Scale (g/N)')
xlabel('Frequency (Hz)')
subplot(2,1,2)
plot(f,phase(Hxx))
ylabel('Phase (rad)')
xlabel('Frequency (Hz)')

figure(4)
subplot(2,1,1)
plot(f,abs(Hxy))
ylabel('Scale (g/N)')
xlabel('Frequency (Hz)')
subplot(2,1,2)
plot(f,phase(Hxy))
ylabel('Phase (rad)')
xlabel('Frequency (Hz)')

figure(5)
subplot(2,1,1)
plot(f,abs(Hxz))
ylabel('Scale (g/N)')
xlabel('Frequency (Hz)')
subplot(2,1,2)
plot(f,phase(Hxz))
ylabel('Phase (rad)')
xlabel('Frequency (Hz)')

figure(6)
subplot(2,1,1)
plot(f,abs(Hyx))
ylabel('Scale (g/N)')
xlabel('Frequency (Hz)')
subplot(2,1,2)
plot(f,phase(Hyx))
ylabel('Phase (rad)')
xlabel('Frequency (Hz)')

figure(7)
subplot(2,1,1)
plot(f,abs(Hyy))
ylabel('Scale (g/N)')
xlabel('Frequency (Hz)')
subplot(2,1,2)
plot(f,phase(Hyy))
ylabel('Phase (rad)')
xlabel('Frequency (Hz)')

figure(8)
subplot(2,1,1)
plot(f,abs(Hyz))
ylabel('Scale (g/N)')
xlabel('Frequency (Hz)')
subplot(2,1,2)
plot(f,phase(Hyz))
ylabel('Phase (rad)')
xlabel('Frequency (Hz)')

figure(9)
subplot(2,1,1)
plot(f,abs(Hzx))
ylabel('Scale (g/N)')
xlabel('Frequency (Hz)')
subplot(2,1,2)
plot(f,phase(Hzx))
ylabel('Phase (rad)')
xlabel('Frequency (Hz)')

figure(10)
subplot(2,1,1)
plot(f,abs(Hzy))
ylabel('Scale (g/N)')
xlabel('Frequency (Hz)')
subplot(2,1,2)
plot(f,phase(Hzy))
ylabel('Phase (rad)')
xlabel('Frequency (Hz)')

figure(11)
subplot(2,1,1)
plot(f,abs(Hzz))
ylabel('Scale (g/N)')
xlabel('Frequency (Hz)')
subplot(2,1,2)
plot(f,phase(Hzz))
ylabel('Phase (rad)')
xlabel('Frequency (Hz)')

% %Emperical TF
% g = etfe(Z);







