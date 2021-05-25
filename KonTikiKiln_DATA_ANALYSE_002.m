clear all
close all
clc

% this program opens a data file and implements some calorimetry
% calcualtions on it to evaluate things like the HRR, mass flows of gases,
% mass loss rate and so on. It then makes a plot of these things and saves
% it. The next step is to make the video aniimation of the TC data.

% USE AT YOUR OWN RISK

% - delay times on the analysers are not considered

% start by chosing a file
[file,path] = uigetfile('*.txt');

%move to the same folder as the file
cd(path)

% then loading it
data = readmatrix(file);

% give the variables some names
timestamp = data(:,1);
test_time  = data(:,2);
v_data = data(:,3:12);
temperatures = data(:,13:17);
test_time  = data(:,18);
mass  = data(:,19);
conc_O2  = data(:,20);
conc_CO  = data(:,21);
conc_CO2  = data(:,22);
DPT_Pa   = data(:,23);
RH_ambient  = data(:,24);
T_cold_trap  = data(:,25);
T_duct  = data(:,26);
T_smoke =  data(:,27);
T_duct2 = data(:,28);
T_ambient   = data(:,29);
Q_OC__O2_CO2_CO  = data(:,30);
HFG_VDC = data(:,31);
TC_Temps = data(:,32:47);

%create a new filename for the outputs
output_filename = [file(1:end-4) '_ANALYSED.mat'];

%%% check to see that ignition (and later flameout) are set

if exist('ignition_row') == 0
    ignition_row = input('When did ignition happen? (vector position)');
    t_ignition = test_time(ignition_row);
else
    t_ignition = test_time(max(ignition_row));
    ignition_row = max(ignition_row);
end

if exist('O2_conc') ==1
    conc_O2 = O2_conc;
    conc_CO = CO_conc;
    conc_CO2 = CO2_conc;
end

save(output_filename);

% create a fighre on which to plot your data
plotfigs = figure('position', [0 0 1500 1000]);


%now change the test time so it is 0 at ignition
test_time = test_time-test_time(max(ignition_row));

% if exist(flameout_row) == 0
%     flameout_row = input('When did flameout happen? (vetcor position)');
%     t_flameout = flameout_row;
%     save(FileName, 'flameout_row', '-append')
% else
%     t_flameout = flameout_row;
% end

%%%% Define a bunch of required variables

%analyser delay times
delay_O2 = 13;
delay_CO2 = 19;
delay_CO2 = 19;

%turn T into units of K
T_ambient = T_ambient+273;
T_duct = T_duct+273;

%MW of species
M_air = 28.8;                                                             % Molecular weight of dry air [kg/kmol]
M_H2O = 18;                                                                 % Molecular weight of H2O [kg/kmol]
M_CO = 28;                                                                  % Molecular weight of CO [kg/kmol]
M_CO2 = 44;                                                                 % Molecular weight of CO2 [kg/kmol]
M_N2 = 28;                                                                  % Molecular weight of N2 [kg/mol]
M_O2 = 32;                                                                  % Molecular weight of O2 [kg/mol]
M_Soot = 12;                                                                % carbon molecular weight [kg/kmol]
M_CH4 = 16;                                                                 % Molecular weight of CH4 [kg/mol]

% universal constants
R = 8.314472;
p_ambient = 101325; %Pa

% energy constants
E_O2 = 13.1;
E_CO2 = 13.3;
E_CO = 17.6;
E_CO_CO2 = 17.63 * 1000; % no idea what this is! assumed to be the energy release from CO->CO2
E_S = 12.3 * 1000;  % 12.3E3[MJ/kg]  Energy release per mass unit of oxygen consumed for combustion of soot to CO2 [kJ/kg]
Delta_H_CO = 283000/28/1000 * 1000; % 2.83E5[kJ/mole] or 10.11 [MJ/kg] Energy release per unit mass of CO consumed in the burning of CO [kJ/mole] Ref.VI
Delta_H_S = 393500/12/1000 * 1000;

% Calorimetry specific constants
alpha = 1.105;

% Apparatus constants
K = 0.86; %probe constant
Dia_Duct = 0.315; %m
Area_Duct = pi * (Dia_Duct^2) / 4;

% Duct flow calculation
rho_air = M_air/1000*101325./(R.*T_duct);
f_Re = 1.08;

% m_Duct_TS = Area_Duct .* K_ ./ f_Re .* ((2 .* p_DPT .* rho_air.^(1/2))); % Mass flow [kg/s] in the test section duct Equation 4 Ref.I
% v_Duct_TS = m_Duct_TS ./ rho_air;                     % Volumetric flow [m3/s] in the test section duct Ref.II

m_duct = 26.54.*Area_Duct.*K./f_Re.*(DPT_Pa./T_duct).^(1/2);
v_duct = m_duct./rho_air.*1000; %l/s

% Mole fractions &c
X_CO_A = conc_CO ./ 10^6;
X_CO2_A = conc_CO2 ./ 10^6;
X_O2_A = conc_O2 ./ 100;

p_H2O_saturation = exp(23.2 - 3816 ./ (-46 + T_ambient));
X0_H2O = RH_ambient .* p_H2O_saturation ./ 100 ./ p_ambient;   % The molar fraction of water vapor Equation (10) Ref.I

X_Soot = 0;

X_H2O = X0_H2O;

X_O2  = X_O2_A  .* (1 - X_H2O - X_Soot);                                         % Molar fraction of oxygen in the exhaust duct from Equation (18) Ref.III
X_CO  = X_CO_A  .* (1 - X_H2O - X_Soot);                                         % Molar fraction of carbon monoxid in the exhaust duct from Equation (18) Ref.III
X_CO2 = X_CO2_A .* (1 - X_H2O - X_Soot);

% averaging the pre-ignition mole fractions
X0_O2_A = mean(X_O2_A(1:ignition_row-1));
X0_CO2_A = mean(X_CO2_A(1:ignition_row-1));
X0_CO_A = mean(X_CO_A(1:ignition_row-1));

% depletion factors and associated stuff

% if O2 and CO2 are measured
phi_O2_CO2 = (X0_O2_A./(1-X_CO2_A) - X_O2_A.*(1-X0_CO2_A))./((1-X_O2_A-X_CO2_A).*X0_O2_A);
% if O2 CO2 and CO are measured
PHI_O2_CO2_CO = (X0_O2_A .* (1 - X_CO2_A - X_CO_A) - ...                    % Equation (20) Ref.I
    X_O2_A .* (1 - X0_CO2_A)) ./ ...
    (X0_O2_A .* (1 - X_O2_A - X_CO2_A - X_CO_A));

m_incoming_air = m_duct ./ (1 + PHI_O2_CO2_CO .* (alpha - 1));             % Mass flow rate of the incomming air [kg/s] Equation (15) Ref.I
M_incoming_air = M_air * (1 - X0_H2O) + M_H2O * X0_H2O;                % Molecular weight of the incomming air [kg/kmol] which should be around 29 Equation (11) Ref.I
m_M__incoming_air_var = m_incoming_air ./ M_incoming_air;
M_air_duct = M_incoming_air;

%% plot gas concentrations to check...
figure(plotfigs)

subplot(2,3,1)
yyaxis left
plot(test_time, conc_O2)
hold on
plot(test_time, conc_CO2./10000)
ylabel('Concentration, %')

yyaxis right
plot(test_time, conc_CO)
xlabel('Time from ignition, s')
ylabel('Concentration, ppm')
legend('[O_2]', '[CO_2]', '[CO]')

%% heat release calcs
% Heat relase by Oxygen Consumption Calorimetry when O2, CO2 and CO are
% measured. From Janssens, M.L. and Parker, W.J. 1992. "Oxygen Consumption
% Calorimetry,"in Heat Release in Fires, edited by V. Babrauskas and S.J. Grayson,Chapter 3: pp. 31-59.

Q_OC__O2_CO2_CO = ((E_O2 .* PHI_O2_CO2_CO - (E_CO - E_O2) .* (1-PHI_O2_CO2_CO) ./ 2 .* X_CO_A ./ X_O2_A) .* ...
    m_incoming_air ./ M_incoming_air .*  M_O2 .* (1 - X0_H2O) * X0_O2_A)*1000;


%%% Calculate the HRR by CDG calorimetry - Note this assumes that all gases
%%% are measured  Brohez et al. (S. Brohez, C. Delvosalle, G. Marlair and
%%% Tewarson, A. 1999. "Soot Generation in Fires: An Important Parameter for Accurate Calculation of Heat Release", Fire Safety Science - Proceedings of the Sixth International Symposium)

Q_CDG__CO2_CO = ((E_CO2 .* M_CO2) .* (m_duct ./ M_incoming_air .* X_CO2_A .* (1 - X0_H2O) - m_incoming_air ./ M_incoming_air .* X0_CO2_A .* (1 - X0_H2O)) + ...
    (E_CO .* M_CO) .* m_duct ./ M_incoming_air .* X_CO_A .* (1 - X0_H2O))*1000;

figure(plotfigs)
subplot(2,3,2)
hold on
box on
plot(test_time, Q_OC__O2_CO2_CO)
plot(test_time, Q_CDG__CO2_CO)
legend('OC O_2, CO, CO_2 (Janssens)', 'CDG CO_2, CO (Brohez et al.)')
xlabel('Time from ignition, s')
ylabel('HRR, kW')

%% mass loss rate

if length(mass)==length(test_time)
    mlr = -diff(mass(1:end))./abs(diff(test_time));
else
    mlr = -diff(mass(1:end-1))./abs(diff(test_time));
end
mlr_smooth = smooth(mlr, 51, 'sgolay', 2)*1000; %g/s

figure(plotfigs)
subplot(2,3,3)
hold on
box on
hold on
plot(test_time(1:end-1), mlr*1000, '.', 'color', [0.7 0.7 0.7])
plot(test_time(1:end-1), mlr_smooth)
% xlim([0 1800])
ylim([-10 10])
legend('Mass loss rate', 'Smoothed mass loss rate')
ylabel('Mass loss rate, g/s')
xlabel('Time from ignition, s')

%% mass flows

m_duct_CO2 = conc_CO2./10^6.*v_duct./1000.*p_ambient*M_CO2./(R.*T_duct); %g/s
m_duct_CO = conc_CO./10^6.*v_duct./1000.*p_ambient*M_CO./(R.*T_duct); %g/s

m_duct_CO2 = m_duct_CO2;
m_duct_CO = m_duct_CO;

% subtract the pre ignition mass

m_duct_CO2 = m_duct_CO2-mean(m_duct_CO2(1:ignition_row-1));
m_duct_CO = m_duct_CO-mean(m_duct_CO(1:ignition_row-1));

figure(plotfigs)
subplot(2,3,4)
hold on
box on
plot(test_time, m_duct_CO2)
plot(test_time, m_duct_CO)
legend('Mass flow CO_2', 'Mass flow CO')
ylabel('Mass flow, g/s')
xlabel('Time from ignition, s')

% yields

yield_CO2 =  m_duct_CO2(1:end-1)./(mlr_smooth(1:end));
yield_CO =  m_duct_CO(1:end-1)./(mlr_smooth(1:end));

figure(plotfigs)
subplot(2,3,5)
hold on
box on
plot(test_time(1:end-1), yield_CO2(1:end))
plot(test_time(1:end-1), yield_CO(1:end))
legend('CO_2', 'CO')
ylabel('Yield, g/g')
xlabel('Time from ignition, s')

% this plots mlr and yelds on the same line
% figure(plotfigs)
% subplot(2,3,6)
% yyaxis left
% plot(test_time(ignition_row:end-1), mlr_smooth(ignition_row:end))
%
% yyaxis right
% hold on
% plot(test_time(ignition_row:end-1), yield_CO2(ignition_row:end))
% plot(test_time(ignition_row:end-1), yield_CO(ignition_row:end))
% legend('CO_2', 'CO')
% ylabel('Yield, g/g')
% xlabel('Time from ignition, s')

%% temperatures

outputvideo = [file(1:end-4) '.mp4'];
output_video = VideoWriter(outputvideo, 'MPEG-4');
output_video.FrameRate = 50;
open(output_video)

figure

for ii = ignition_row:length(temperatures)

    figure('visible','off')

    plot(TC_Temps(ii,1:4), [1 2 3 4])
    hold on
    plot(TC_Temps(ii,5:8), [1 2 3 4])
    plot(TC_Temps(ii,9:12), [1 2 3 4])
    plot(TC_Temps(ii,13:16), [1 2 3 4])
    
    xlim([0 1000])
    xlabel('Temperature, °C')
    ylabel('Height, mm')
    legend('Position 1', 'Position 2', 'Position 3', 'Position 4')
    
    % subplot(2,2,1)
    % plot(TC_Temps(ii,1:4), [1 2 3 4])
    % xlim([0 1000])
    % xlabel('Temperature, °C')
    % ylabel('Height, mm')
    % subplot(2,2,2)
    % plot(TC_Temps(ii,5:8), [1 2 3 4])
    % xlim([0 1000])
    % xlabel('Temperature, °C')
    % ylabel('Height, mm')
    % subplot(2,2,3)
    % plot(TC_Temps(ii,9:12), [1 2 3 4])
    % xlim([0 1000])
    % xlabel('Temperature, °C')
    % ylabel('Height, mm')
    % subplot(2,2,4)
    % plot(TC_Temps(ii,13:16), [1 2 3 4])
    % xlim([0 1000])
    % xlabel('Temperature, °C')
    % ylabel('Height, mm')
    
    frame = getframe(gcf);
    
    writeVideo(output_video, frame);
    close
end

close(output_video)

%% saving

save(output_filename, '-regexp', '^(?!(file|exptname|jj|output_video|plotfigs)$).')

saveas(gcf, [file(1:end-4) '.png'])


