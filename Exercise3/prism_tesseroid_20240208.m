% potential of a prism
%
% status: 2022-11-14
%
clear all; close all; clc; format long g;

% general constants
G       = 6.6743*10^(-11);     % Newtons gravitational constant
density = 2670;                % Density of the shell
R = 6371000;

% size of the prism
% lower, left, front
xc(1) = -500
yc(1) = -500
zc(1) = -500
% uper, right,top
xc(2) = -xc(1)
yc(2) = -yc(1)
zc(2) = -zc(1)

% computation popint P
zP = [zc(1)-10000:1:zc(2)+10000];
nP = length(zP);
xP = ones(1,nP)*(xc(2)+xc(1))/2;
yP = ones(1,nP)*(yc(2)+yc(1))/2;





% Tesseroid
% Dimension of the Tesseroids
dr = zc(2)-zc(1)
dphi = (xc(2)-xc(1))/R
dlam = (yc(2)-yc(1))/R
% Square of the dimensions (Tesseroidformula)
dr2 = dr*dr;
dphi2 = dphi*dphi;
dlam2 = dlam*dlam;

% Coordinates of the geometrical center of the Tesseroid
r_0  = R;
phi_0 = 0;
lam_0 = 0;
sin_phi_0 = sin(phi_0);
cos_phi_0 = cos(phi_0);
sin_lam_0 = sin(lam_0);
cos_lam_0 = cos(lam_0);

r_02 = r_0*r_0;
sin_phi_02 = sin_phi_0*sin_phi_0;
konst1 = G*density*dr*dphi*dlam;

phi_P = 0;
lam_P = 0;
sin_phi_P = sin(phi_P);
cos_phi_P = cos(phi_P);
sin_lam_P = sin(lam_P);
cos_lam_P = cos(lam_P);

for n = 1:nP
    u(n) = 0;
    V(n) = 0;
    punktmasse(n)=0;
for i = 1:2
for j = 1:2
for k = 1:2
   vz = (-1)^(i+j+k);
   dx = xP(n)-xc(i);
   dy = yP(n)-yc(j);
   dz = zP(n)-zc(k);
   Wijk = sqrt( dx^2 + dy^2 + dz^2  );
   u(n) = u(n) - vz*(dx*dy*abs(log(dz+Wijk)) + dy*dz*abs(log(dx+Wijk)) + dx*dz*abs(log(dy+Wijk)) );
   if abs(dx) > 10^(-10)
       u(n) = u(n) + 0.5 * vz * dx^2 * atan(dy*dz/(dx*Wijk));
   end
   if abs(dy) > 10^(-10)
       u(n) = u(n) + 0.5 * vz *  dy^2 * atan(dz*dx/(dy*Wijk));
   end
   if abs(dz) > 10^(-10)
       u(n) = u(n) + 0.5 * vz *  dz^2 * atan(dy*dx/(dz*Wijk));
   end
end
end
end
u(n) = u(n)*G*density;


rp   = R + zP(n);
  rp2 = rp*rp;
  rPr0 = rp*r_0;
  rP2_r02 = rp2 + r_02;
  sin_dlam = sin_lam_0*cos_lam_P - cos_lam_0*sin_lam_P;
  cos_dlam = cos_lam_0*cos_lam_P + sin_lam_0*sin_lam_P;   

  % see Heck and Seitz (2007), equation (26) sphereical distance betwen P and Q
  cos_psi  = sin_phi_0*sin_phi_P + cos_phi_0*cos_phi_P*cos_dlam;
  sin_psi2 = (1-cos_psi*cos_psi);

  spPsp0 = sin_phi_P*sin_phi_0;
  cpPcp0 = cos_phi_P*cos_phi_0;

  % see Heck and Seitz (2007), equation (25) Euclidean distance betwen P and Q
  ell2 = (rP2_r02 - 2*rPr0*cos_psi);
  ell  = sqrt(ell2)
  ell5 = ell2*ell2*ell;

  % see Heck and Seitz (2007), equation (25) K_000
  K000 = r_02*cos_phi_0/ell
  
  % see Heck and Seitz (2007), equation (43) K_200
  K200 = (rp2*cos_phi_0/ell5)*(2*ell2 - 3*r_02*sin_psi2);

  % see Heck and Seitz (2007), equation (44) K_020
  K020 = (r_02/ell5)*(-cos_phi_0*rP2_r02*(rP2_r02 -...
          rPr0*spPsp0)+(rp2*r_02*cos_phi_0)*...
          (sin_phi_P*sin_phi_P*(3-sin_phi_02)-(cos_phi_P*cos_phi_P)*...
          (2-sin_phi_02)*cos_dlam*cos_dlam)+(rPr0*cos_phi_P)*...
          (3-sin_phi_02)*(rP2_r02 - 2*rPr0*spPsp0)*cos_dlam);
    
  % see Heck and Seitz (2007), equation (45) K_002
  K002 = -(rPr0*r_02*cpPcp0*cos_phi_0/ell5)*...
         (ell2*cos_dlam - 3*rPr0*cpPcp0*sin_dlam*sin_dlam);
   
  % see Heck and Seitz (2007), equation (24) Gravitationspotential V
  V(n) = konst1*(K000 + (K200*dr2 + K020*dphi2 + K002*dlam2)/24);

  punktmasse(n) = konst1*K000;

diff(n) = abs(u(n)-V(n));
end


figure
semilogy(zP,u, zP,punktmasse, zP,V,  'LineWidth',1)
legend( 'prism', 'pointmass', 'tesseroid' )

