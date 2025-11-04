% uniaxial tension - perpendicular
clear
clc
close all

tic
% setting up loading
edot = 0.01;       % loading rate. unit: 1/s
emax = 2.0;         % maximal stretch.
dt = 0.01;             % time increment used for simulation. unit: s
tend = emax/edot;
Nstep2 = tend/dt;

% save loading path information
load_path.dt = dt;
load_path.edot = edot;
load_path.emax = emax;

% load initial state
load initial_state.mat

% initialize temp variables
Resi = zeros(7,1);
Reso = 0.1*ones(3,1);
Ki = zeros(7,7);
KF = zeros(3,3);
KQ = zeros(3,7);
GF = zeros(7,3);
countermax = 100;
hitmax = false;

Fp = F_eq;
Fvp = Fv_eq;
Fep = Fe_eq;
qp = q_eq;
Sp = S_eq;
Tn = T_eq;   

Fn = Fp;
Fvn = Fvp;
Fen = Fep;
qn = qp;
Sn = Sp;

% initialize ouput variables
time_save = [];         % time
Fref22_save = [];       % F_22 component of deformation gradient, scaled by Feq
sref22_save = [];       % sigma_22 component of Cauchy stress
qref_save = [];         % director (in radian)
Sref_save = [];         % order parameter
Fvref22_save = [];      % Fv_22 component of viscous strain
F22_save = [];          
Fref12_save = [];       % F_12 component of deformation gradient, scaled by Feq
F12_save = [];

% time marching
for istep = 1:2*Nstep2+1
    % update time and loading
    tn = istep*dt;
    if istep < Nstep2
        Fn(2,2) = Fn(2,2) + edot*dt;
    elseif sigman(2,2) >= 0
        Fn(2,2) = Fn(2,2) - edot*dt;
    else
        hitmax = true;
    end

    Reso = 0.1*ones(2,1);
    counter1 = 0;
    while (sqrt(Reso.'*Reso) > 1e-11) && ~hitmax
        counter1 = counter1 + 1;
        [Rd,RS,RFv11,RFv12,RFv21,RFv22,RFv33] = residual_MS(mp,dt,Tn,Fn,Fen,Fvn,qn,Sn,Fp,Fvp,qp,Sp);
        Resi(1) = Rd;
        Resi(2) = RS;
        Resi(3) = RFv11;
        Resi(4) = RFv22;
        Resi(5) = RFv33;
        Resi(6) = RFv12;
        Resi(7) = RFv21;

        counter2 = 0;
        if Sn < 1e-6
            Sn = 1e-6;
        end
        while (sqrt(Resi.'*Resi) > 1e-13) && ~hitmax
            counter2 = counter2 + 1;
            [Ki,Ko] = tangents_MS(mp,dt,Tn,Fn,Fen,Fvn,qn,Sn,Fp,Fep,Fvp,qp,Sp);
            delta = - inv(Ki)*Resi;
            qn = qn + delta(1);
            Sn = Sn + delta(2);
            Fvn(1,1) = Fvn(1,1) + delta(3);
            Fvn(2,2) = Fvn(2,2) + delta(4);
            Fvn(3,3) = Fvn(3,3) + delta(5);
            Fvn(1,2) = Fvn(1,2) + delta(6);
            Fvn(2,1) = Fvn(2,1) + delta(7);
            Fen = Fn*inv(Fvn);
            sigman = stress_MS(mp,Fn,Fen,qn,Sn);
            [Rd,RS,RFv11,RFv12,RFv21,RFv22,RFv33] = residual_MS(mp,dt,Tn,Fn,Fen,Fvn,qn,Sn,Fp,Fvp,qp,Sp);
            Resi(1) = Rd;
            Resi(2) = RS;
            Resi(3) = RFv11;
            Resi(4) = RFv22;
            Resi(5) = RFv33;
            Resi(6) = RFv12;
            Resi(7) = RFv21;
            if counter2 == countermax
                fprintf('%f, %f, %f, %f, %f, %f, %f\n',Resi(1),Resi(2),Resi(3),Resi(4),Resi(5),Resi(6),Resi(7));
                hitmax = true;
            end
        end

        % outer loop
        Reso(1) = sigman(1,1);
        Reso(2) = sigman(1,2);
        [Ki,Ko] = tangents_MS(mp,dt,Tn,Fn,Fen,Fvn,qn,Sn,Fp,Fep,Fvp,qp,Sp);
        Ko_new = [Ko(1,1) Ko(1,3);Ko(3,1) Ko(3,3)];
        deltao = - inv(Ko_new)*Reso;
        Fn(1,1) = Fn(1,1) + deltao(1);
        Fn(1,2) = Fn(1,2) + deltao(2);
        Fn(3,3) = 1/(Fn(1,1)*Fn(2,2));
        Fen = Fn*inv(Fvn);
        sigman = stress_MS(mp,Fn,Fen,qn,Sn);
        Reso(1) = sigman(1,1);
        Reso(2) = sigman(1,2);

        if counter1 == countermax
            fprintf('%f, %f, %f\n',Reso(1),Reso(2));
            hitmax = true;
        end
    end

    % calculate driving force - post-processing
    Fn_dot = (Fn - Fp)/dt;
    Ln = Fn_dot*inv(Fn);
    Wn = 0.5*(Ln-Ln.');

    % update 'past' variables
    Fp = Fn;
    Fvp = Fvn;
    qp = qn;
    Sp = Sn;

    % save for output
    F22_save = [F22_save;Fn(2,2)];
    Fref = inv(F_eq)*Fn;
    time_save = [time_save;tn];
    Fref22_save = [Fref22_save;Fref(2,2)];
    sref22_save = [sref22_save;sigman(2,2)];
    qref_save = [qref_save;qn];
    Sref_save = [Sref_save;Sn];
    Fvref22_save = [Fvref22_save;Fvn(2,2)];
    F12_save = [F12_save;Fn(1,2)];
    Fref12_save = [Fref12_save;Fref(1,2)];

end

figure()
subplot(2,2,1)
plot(Fref22_save,sref22_save,'r-','linewidth',2);
grid on
xlabel('F_{22}','fontsize',15)
ylabel('\sigma_{22}','fontsize',15)
title('\sigma_{22}-F_{22}')

subplot(2,2,2)
plot(Fref22_save,Fvref22_save,'r-','linewidth',2);
grid on
xlabel('F_{22}','fontsize',15)
ylabel('F^v_{22}','fontsize',15)
title('F^v_{22}-F_{22}')

subplot(2,2,3)
plot(Fref22_save,qref_save,'r-','linewidth',2);
grid on
xlabel('F_{22}','fontsize',15)
ylabel('\theta','fontsize',15)
title('\theta-F_{22}')

subplot(2,2,4)
plot(Fref22_save,Sref_save,'r-','linewidth',2);
grid on
xlabel('F_{22}','fontsize',15)
ylabel('S','fontsize',15)
title('S-F_{22}')

clearvars -except time_save Fref22_save sref22_save qref_save Sref_save...
    Fvref22_save Fref12_save mp load_path
save result_uniaxial_tension_perp time_save Fref22_save sref22_save...
    qref_save Sref_save Fvref22_save Fref12_save mp load_path
toc