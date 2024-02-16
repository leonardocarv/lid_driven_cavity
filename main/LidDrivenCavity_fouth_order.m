clear all; close all; clc

% Model parameters

u_wall = 1;    % sliding wall velocity (m/s)
Re     = 100;  % Reynolds Number
Pr     = 1;    % Prandtl Number
Sc     = 1;    % Schmidt Number
Ri     = 0.01; % Richardson Number
N      = 100;  % Buoyant Ratio Number

np     = 81;   % Number of grid elements in x and y directions


% Initiating model specific parameters

omega_err      = 1;
omega_err_past = 1;
theta_err      = 1;
psi_err        = 1;
err            = 1e-6;
beta           = 1.7;

x = linspace(0,1,np);
y = linspace(0,1,np);

[X, Y] = meshgrid(x,y);

omega = zeros(length(y),length(x));   %vorticity
psi   = zeros(length(y),length(x));   %streamfunction
theta = zeros(length(y),length(x));   %temperature
phi   = zeros(length(y),length(x));   %mass


clear x y c i dx np

it = 0;
it2 = 0;

dx = X(1,2) - X(1,1);

ds = 1/((1/(dx*dx)) + (1/(dx*dx)));

dt_Re = 0.5*ds*Re;

dt = dt_Re/1;

while omega_err > err
  tic;
  psi_err = 1;
  it = it + 1;

  dy1 = (Y(2,1)-Y(1,1));
  dye = (Y(end,1)-Y(end-1,1));

  % upper boundary

  omega(end,:) = (-2/(dy1*dy1))*(psi(end-1,:)) - u_wall*2/dy1;
  theta(end,:) = 1;
  phi(end,:)   = 1;

  %bottom boundary

  omega(1,:) = (-2/(dye*dye))*(psi(2,:));
  theta(1,:) = 0;
  phi(1,:)   = 0;


  %left boundary

  omega(:,1) = -2*(psi(:,2))/(dy1*dy1);
  theta(:,1) = theta(:,2);
  phi(:,1)   = phi(:,2);

  %right boundary

  omega(:,end) = -2*(psi(:,end-1))/(dy1*dy1);
  theta(:,end) = theta(:,end-1);
  phi(:,end)   = phi(:,end-1);

  omega_aux = omega;

  L = zeros(length(X),length(Y),5);
  L_phi = zeros(length(X),length(Y),5);
  L_theta = zeros(length(X),length(Y),5);

  omg_aux = zeros(length(X),length(Y),5);
  omg_aux(:,:,1) = omega;
  phi_aux = zeros(length(X),length(Y),5);
  phi_aux(:,:,1) = phi;
  theta_aux = zeros(length(X),length(Y),5);
  theta_aux(:,:,1) = theta;

  for k = 1:5
    for i = 2:length(Y)-1
      for j = 2:length(X)-1

        dx = (X(i,j+1) - X(i,j));
        dy = (Y(i+1,j) - Y(i,j));

        %advection of vorticity, temperature and concentration

        [omg_advx,omg_advy] = doAdvection(psi(i+1,j),psi(i-1,j),psi(i,j+1),psi(i,j-1),...
        omega(i,j+1),omega(i,j-1),omega(i+1,j),omega(i-1,j),dx,dy); % vorticity

        [phi_advx,phi_advy] = doAdvection(psi(i+1,j),psi(i-1,j),psi(i,j+1),psi(i,j-1),...
        phi(i,j+1),phi(i,j-1),phi(i+1,j),phi(i-1,j),dx,dy); % concentration

        [theta_advx,theta_advy] = doAdvection(psi(i+1,j),psi(i-1,j),psi(i,j+1),psi(i,j-1),...
        theta(i,j+1),theta(i,j-1),theta(i+1,j),theta(i-1,j),dx,dy); % temperature

        % diffusion of vorticity, temperature and concentration
        [omg_diffx,omg_diffy] = doDiffusion(omega(i,j),omega(i,j+1),omega(i,j-1),...
        omega(i+1,j),omega(i-1,j),dx,dy); % vorticity

        [phi_diffx,phi_diffy] = doDiffusion(phi(i,j),phi(i,j+1),phi(i,j-1),...
        phi(i+1,j),phi(i-1,j),dx,dy); % concentration

        [theta_diffx,theta_diffy] = doDiffusion(theta(i,j),theta(i,j+1),theta(i,j-1),...
        theta(i+1,j),theta(i-1,j),dx,dy); % temperature

        % time evolution
        L(i,j,k) = (-1*omg_advx + omg_advy - (1/Re)*(-1*omg_diffx - omg_diffy)...
        + (Ri*Ri)*(((theta(i,j+1) - theta(i,j-1))/(2*dx)) + N*((psi(i,j+1) - psi(i,j-1))/(2*dy)))); %vorticity

        L_phi(i,j,k) = ((-1*phi_advx + phi_advy + (1/(Re*Sc))*(phi_diffx + phi_diffy))); % concentration

        L_theta(i,j,k) = ((-1*theta_advx + theta_advy + (1/(Re*Pr))*(theta_diffx + theta_diffy))); % temperature

        clear omg_advx omg_advy omg_diffx omg_diffy theta_advx theta_advy theta_diffx theta_diffy ...
        phi_advx phi_advy phi_diffx phi_diffy

      end
    end

    omg_aux(:,:,k+1)   = doTimeEvolution(L,omg_aux,k,dt);
    phi_aux(:,:,k+1)   = doTimeEvolution(L_phi,phi_aux,k,dt);
    theta_aux(:,:,k+1) = doTimeEvolution(L_theta,theta_aux,k,dt);
    omega              = omg_aux(:,:,k+1);
    phi                = phi_aux(:,:,k+1);
    theta              = theta_aux(:,:,k+1);

  end

  while psi_err > err
    it2 = it2 + 1;
    psi_aux = psi;
    for i = 2:length(Y)-1
      for j = 2:length(X)-1
        dx = abs(X(i,j+1)-X(i,j));
        dy = abs(Y(i+1,j)-Y(i,j));
        b = 1/(dx*dx);
        c = 1/(dy*dy);
        a = b + c;

        A = b*(psi(i+1,j) + psi(i-1,j)) + c*(psi(i,j+1)+psi(i,j-1));
        B = (1/12)*(omega(i+1,j)+omega(i,j+1)+8*omega(i,j)+omega(i-1,j)+omega(i,j-1));
        C = (1/12)*a*(psi(i+1,j+1) - 2*psi(i,j+1)+psi(i-1,j+1)-2*psi(i+1,j)...
        -2*psi(i-1,j)+psi(i+1,j-1)-2*psi(i,j-1)+psi(i-1,j-1));
        psi(i,j) = (3/5)*(beta/a)*(A+B+C)+(1-beta)*psi(i,j);
      end
    end

    psi_err = abs(max(abs(psi_aux(:))) - max(abs(psi(:))));
    clear dx dy a b c A B C
  end

  elapsed_time = toc;
  omega_err = abs(max(max(abs(omega) - abs(omega_aux))));

  fprintf('%-6g %E %10.7f %10.7f\n',it,omega_err,dt,elapsed_time)

end


