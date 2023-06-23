%% TRY TO RETRIEVE THE "IMAGE" TRANSFORM BY INVERSE FOURIER TRANSFORM

%% 1D SIGNAL
% f(x) = g(phi(x))
% f(x) = 1/nX * int_Y{G(Y) exp(i.Y.phi(x)) dY}
% F(X) = int_x{ 1/nX * int_Y{G(Y) exp(i[Y.phi(x)-X.x]) dY} dx} = int_Y{G(Y) P(X,Y) dY}
% with P(X,Y) = 1/nX * int_x{exp(i[Y.phi(x)-X.x]) dx}
% P(X,Y) = 1/nX * int_x{Z(x).^Y .* exp(-i.X.x) dx} = 1/nX * int_x{V(Y,x) .* exp(-i.X.x) dx}
% with Z(x) = exp(1i*phi(x)) 
% and V(Y,x) = Z(x).^Y is an Vandermonde matrix

nX = 100 ;
x = 0:nX-1 ;
d = 2*pi/(nX) ;
X = x ; 

g = exp(-(9*(x-nX/2)/nX).^2) ;
u = ...1*nX*((14/100)*(2*(x/nX-.5)).^2 + 15/1000) ...
    1/100*nX*sin(2*pi*(x/nX-.5)*21) ...
    ;
phi = x + u ;
f = interp1([x-nX x x+nX],[g g g],phi,'cubic') ;

F = fft(f) ;
G = fft(g) ;
F = fftshift(F) ; G = fftshift(G) ; X = [X(X>=nX/2)-nX , X(X<nX/2)] ; % fftshift

V = exp(1i*d.*phi).^X(:) ;
fb = 1/nX * G*V ;
err_f = norm(fb-f)/norm(f)

E = exp(-1i*d*X.*x(:)) ;
P = 1/nX*V*E ;
Fa = (G*P) ;
err_F = norm(Fa-F)/norm(F)

clf ; plot(x,phi-x) ;
clf ; plot(x,f) ; plot(x,g) ; plot(x,real(fb),':') ;
%clf ; plot3(X,real(F),imag(F)) ; plot3(X,real(G),imag(G)) ; plot3(X,real(Fa),imag(Fa),':') ;
%clf ; imagesc(real(P))
%%
Ga = 1/(nX*1+0) * (fft(g)) ; y = x ;
Ga = fftshift(Ga) ; y = X ;
z = zeros(nX-1,nX) ;
for mm=1:nX
    h = exp(1i*d*(mm-1)).^y ;
    p = Ga.*h ;
    p(y==0) = p(y==0)-f(mm) ;
    z(:,mm) = roots(flip(p)) ;
end

mm = 19 ; nX/2 ; clf ; axis equal ; plot3(real(z(:,mm)),imag(z(:,mm)),ones(nX-1,1).*x(mm)/nX,'.')
%clf ; axis equal ;
zt = exp(1i*d*phi) ;
plot(real(zt),imag(zt),'+')
plot(real(zt(mm)),imag(zt(mm)),'o')
plot(real(exp(1i*d*u(mm))),imag(exp(1i*d*u(mm))),'o')

