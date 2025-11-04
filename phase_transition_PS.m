% test and phase transition simulation

clear
clc
close all

tic
% material properties
mp.etaD0 = 40;      % director viscosity. unit: MPa*s
mp.etaS = 40;       % order viscosity. unit: MPa*s
mp.etaN = 40;       % bilk network viscosity. unit: MPa*s
mp.ks = 0.7;        % nonlinear coefficient in the director viscosity
mp.a0 = 0.3;        % parameter of nematic energy. unit: MPa/K
mp.a1 = 1e-3;      % parameter of nematic energy. unit: MPa
mp.a2 = 40.0;         % parameter of nematic energy. unit: MPa
mp.Tni = 373;       % transition temperature. unit: K
mp.mueq = 1.0;      % equilibrium modulus. unit: MPa
mp.muneq = 1.0;     % non-equilibrium modulus. unit: MPa

% specify thermal path
dt = 1.0;             % time increment used for simulation. unit: s
Tdot = 1/60;        % cooling rate. unit: K/s
Thigh = 383;        % iniital temperature. unit: K
Tlow = 323;         % ultimate temperature. unit:K

Tanneal = 300;      % annealing time before cooling. unit: K
Tstand = 1;     % annealing time after cooling. unit: K
Tload = (Thigh-Tlow)/Tdot;
Nstep = Tload/dt + Tanneal/dt + Tstand/dt;

% save thermal path information
thermal_path.dt = dt;
thermal_path.Tdot = Tdot;
thermal_path.Thigh = Thigh;
thermal_path.Tlow = Tlow;
thermal_path.Tanneal = Tanneal;
thermal_path.Tstand = Tstand;

% initial condition
Fp = eye(3);
Fvp = Fp;
Fep = Fp*inv(Fvp);
qp = 0.5/180*pi;
Sp = 1e-5;

Fn = Fp;
Fvn = [1 0 0;0 1+1e-6 0;0 0 1];
Fen = Fn*inv(Fvn);
qn = qp;
Sn = Sp;

Resi = zeros(7,1);
Reso = 0.1*ones(3,1);
Ki = zeros(7,7);
KF = zeros(3,3);
KQ = zeros(3,7);
GF = zeros(7,3);
countermax = 40;
hitmax = false;

% initialize output variables
time_save = [];     % time
temp_save = [];     % temperature
F11_save = [];      % F_11 component of deformation gradient
F22_save = [];      % F_22 component of deformation gradient
F33_save = [];
S_save = [];        % order parameter
q_save = [];        % director (in radian)
Fv11_save = [];
Fv22_save = [];


% time marching
for istep = 1:Nstep
    % update time and temperature
    tn = istep*dt;
    if tn < Tanneal
        Tn = Thigh;
    elseif (tn >= Tanneal) && (tn < Tanneal+Tload)
        Tn = Tn - Tdot*dt;
    elseif ( tn >= Tanneal+Tload) && (tn <= Tanneal+Tload+Tstand)
        Tn = Tn;
    end

    Reso = 1*ones(3,1);
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
        if Sn < 1e-10
            Sn = 1e-10;
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
                hitmax = True;
            end
        end
        % if Sn < 1e-10
        %     Sn = 1e-10;
        % end

        % outer loop
        Reso(1) = sigman(1,1);
        Reso(2) = sigman(2,2);
        Reso(3) = sigman(1,2);
        [Ki,Ko] = tangents_MS(mp,dt,Tn,Fn,Fen,Fvn,qn,Sn,Fp,Fep,Fvp,qp,Sp);
        deltao = - inv(Ko)*Reso;
        Fn(1,1) = Fn(1,1) + deltao(1);
        Fn(2,2) = Fn(2,2) + deltao(2);
        Fn(1,2) = Fn(1,2) + deltao(3);
        Fn(3,3) = 1/(Fn(1,1)*Fn(2,2));
        Fen = Fn*inv(Fvn);
        sigman = stress_MS(mp,Fn,Fen,qn,Sn);
        Reso(1) = sigman(1,1);
        Reso(2) = sigman(2,2);
        Reso(3) = sigman(1,2);

        if counter1 == countermax
            fprintf('%f, %f, %f\n',Reso(1),Reso(2),Reso(3));
            hitmax = true;
        end
    end
    % update 'past' variables
    Fp = Fn;
    Fvp = Fvn;
    qp = qn;
    Sp = Sn;

    % save for output
    time_save = [time_save;tn];
    temp_save = [temp_save;Tn];
    S_save = [S_save;Sn];
    q_save = [q_save;qn];
    F11_save = [F11_save;Fn(1,1)];
    F22_save = [F22_save;Fn(2,2)];
    F33_save = [F33_save;Fn(3,3)];
    Fv11_save = [Fv11_save;Fvn(1,1)];
    Fv22_save = [Fv22_save;Fvn(2,2)];
end
% obtain equilibrium state at the time end
F_eq = Fn;
Fv_eq = Fvn;
Fe_eq = Fn/Fv_eq;
S_eq = Sn;
T_eq = Tn;
q_eq = qn;

figure()
subplot(2,3,1)
plot(time_save,temp_save,'k-','linewidth',2);
grid on
xlabel('time t(s)','fontsize',15)
ylabel('temperature T(K)','fontsize',15)
title('T-t')

subplot(2,3,2)
plot(temp_save,S_save,'r-','linewidth',2)
grid on
xlabel('temperature T(K)','fontsize',15)
ylabel('order S','fontsize',15)
title('S-T')

subplot(2,3,3)
plot(temp_save,q_save,'r-','linewidth',2)
grid on
xlabel('temperature T(K)','fontsize',15)
ylabel('director \theta','fontsize',15)
title('\theta-T')

subplot(2,3,4)
plot(temp_save,F11_save,'r-','linewidth',2)
grid on
xlabel('temperature T(K)','fontsize',15)
ylabel('F_{11}','fontsize',15)
title('F_{11}-T')

subplot(2,3,5)
plot(temp_save,F22_save,'r-','linewidth',2)
grid on
xlabel('temperature T(K)','fontsize',15)
ylabel('F_{22}','fontsize',15)
title('F_{22}-T')

subplot(2,3,6)
plot(temp_save,Fv11_save,'r-','linewidth',2)
grid on
xlabel('temperature T(K)','fontsize',15)
ylabel('F^v_{11}','fontsize',15)
title('F^v_{11}-T')

figure()
subplot(2,3,1)
plot(time_save,temp_save,'k-','linewidth',2);
grid on
xlabel('time t(s)','fontsize',15)
ylabel('temperature T(K)','fontsize',15)
title('T-t')

subplot(2,3,2)
plot(time_save,S_save,'b-','linewidth',2)
grid on
xlabel('time t(s)','fontsize',15)
ylabel('order S','fontsize',15)
title('S-T')

subplot(2,3,3)
plot(time_save,q_save,'b-','linewidth',2)
grid on
xlabel('time t(s)','fontsize',15)
ylabel('director \theta','fontsize',15)
title('\theta-T')

subplot(2,3,4)
plot(time_save,F11_save,'b-','linewidth',2)
grid on
xlabel('time t(s)','fontsize',15)
ylabel('F_{11}','fontsize',15)
title('F_{11}-T')

subplot(2,3,5)
plot(time_save,F22_save,'b-','linewidth',2)
grid on
xlabel('time t(s)','fontsize',15)
ylabel('F_{22}','fontsize',15)
title('F_{22}-T')

subplot(2,3,6)
plot(time_save,Fv11_save,'b-','linewidth',2)
grid on
xlabel('time t(s)','fontsize',15)
ylabel('F^v_{11}','fontsize',15)
title('F^v_{11}-T')

% load pure_lc.mat
% figure()
% plot(Expression1(:,1)-273,Expression1(:,2));hold on
% plot(temp_save-273,S_save)

clearvars -except mp thermal_path time_save temp_save F11_save F22_save...
    S_save q_save Fv11_save Fv22_save F_eq Fe_eq Fv_eq S_eq q_eq T_eq

save 'initial_state.mat' mp thermal_path F_eq Fe_eq Fv_eq S_eq q_eq T_eq
save 'results_phase_transition' mp thermal_path time_save temp_save F11_save...
    F22_save S_save q_save Fv11_save Fv22_save
toc