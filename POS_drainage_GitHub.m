% Post Oak Savanna Drainage
% By mingxiu Wang on 06/28/24

%% unit gradient approach for site 1 (savanna grassland) using sensor MP data
daily_mean_MP_site1_cm = readtimetable("daily_mean_MP_site1_cm.xlsx");
% VG parameters
VGparameter = readmatrix('SoilTexture.xlsx', 'Sheet', 'Sheet2', 'Range', 'C3:G22');

% For site1
theta_r = VGparameter(1:5,1);
theta_s = VGparameter(1:5,2);
alpha = VGparameter(1:5,3);
n = VGparameter(1:5,4);
Ks = VGparameter(1:5,5); % cm/day
m_VG = 1-1./n;
% L=-1;
L=0.5;

% Calculate K
K_psi_store = zeros(height(daily_mean_MP_site1_cm),5);

% try K(Se), only when L= 0.5, same to the K(h)!!!!!!!!!!!!!!
for i = 1:5
    alpha_psi = abs(alpha(i).*table2array(daily_mean_MP_site1_cm(:,i))); % to solve the complex double results, should be absolute value!
%     S = (1/(1+alpha_psi.^n(i))).^m_VG(i);
    S = (1+(alpha_psi).^n(i)).^(-m_VG(i));
    K_psi = Ks(i).*S.^L.*(1-(1-S.^(1./m_VG(i))).^m_VG(i)).^2;
    % Store the K values
    K_psi_store(:,i) = K_psi; % cm/day
end

q_UG_site1 = -K_psi_store(:,5); % K at 100 cm (q=-K()) -- cm/day

%% Arithmetic mean K for site 1 (savanna grassland) using sensor MP data
K_mean_store = zeros(height(daily_mean_MP_site1_cm),4);
for i = 1:4
    K_mean = mean(K_psi_store(:,i:5),2);
    K_mean_store(:,i) = K_mean;
end

h_gradient_20_100 = (daily_mean_MP_site1_cm.corrected_MP_site1_cm5 - daily_mean_MP_site1_cm.corrected_MP_site1_cm1)/(-80);
h_gradient_40_100 = (daily_mean_MP_site1_cm.corrected_MP_site1_cm5 - daily_mean_MP_site1_cm.corrected_MP_site1_cm2)/(-60);
h_gradient_60_100 = (daily_mean_MP_site1_cm.corrected_MP_site1_cm5 - daily_mean_MP_site1_cm.corrected_MP_site1_cm3)/(-40);
h_gradient_80_100 = (daily_mean_MP_site1_cm.corrected_MP_site1_cm5 - daily_mean_MP_site1_cm.corrected_MP_site1_cm4)/(-20);
h_gradient = [h_gradient_20_100 h_gradient_40_100 h_gradient_60_100 h_gradient_80_100];

q_method2_site1 = -K_mean_store.*(1+h_gradient);

%% Layer-based K with harmonic mean for site 1 (savanna grassland) using sensor MP data
layerK_store = zeros(height(daily_mean_MP_site1_cm),4);
% L = 20; % Length of each layer is 20 cm
for i = 1:4
    layerK = mean(K_psi_store(:,[i,i+1]),2);
    layerK_store(:,i) = layerK;
end

K_eff_store = zeros(height(daily_mean_MP_site1_cm),3);
% Compute harmonic mean for depths 20-100 cm, 40-100 cm, and 60-100 cm
for i = 1:3
    
    K_eff = harmmean(layerK_store(:,i:4),2);
    K_eff_store(:,i) = K_eff;
end

q_method3_site1 = -K_eff_store.*(1+h_gradient(:,1:3));%cm/day

% Create a new timetable combining the time and q_UG_site4
q_combined_site1 = timetable(daily_mean_MP_site1_cm.Time, q_UG_site1,q_method2_site1,q_method3_site1);
% Rename the variable for clarity (optional)
q_combined_site1.Properties.VariableNames = {'q_UG','q_method2','q_method3'};

sumQ1_site1 = sum(q_UG_site1);
sumQ2_site1 = sum(q_method2_site1);
sumQ3_site1 = sum(q_method3_site1);

annual_Q1_site1 = sumQ1_site1/(height(q_UG_site1)/365);
annual_Q2_site1 = sumQ2_site1/(height(q_method2_site1)/365);
annual_Q3_site1 = sumQ3_site1/(height(q_method3_site1)/365);
%% K with harmonic mean for site 1 (savanna grassland) using sensor MP data
K_eff_store = zeros(height(daily_mean_MP_site1_cm),4);
% Compute harmonic mean for depths 20-100 cm, 40-100 cm, and 60-100 cm
for i = 1:4
    K_eff = harmmean(K_psi_store(:,i:5),2);
    K_eff_store(:,i) = K_eff;
end

q_method4_site1 = -K_eff_store.*(1+h_gradient);%cm/day

% Create a new timetable combining the time and q_UG_site4
q_combined_site1 = timetable(daily_mean_MP_site1_cm.Time, q_UG_site1,q_method2_site1,q_method3_site1,q_method4_site1);
% Rename the variable for clarity (optional)
q_combined_site1.Properties.VariableNames = {'q_UG','q_method2','q_method3','q_method4'};

sumQ1_site1 = sum(q_UG_site1);
sumQ2_site1 = sum(q_method2_site1);
sumQ3_site1 = sum(q_method3_site1);
sumQ4_site1 = sum(q_method4_site1);

annual_Q1_site1 = sumQ1_site1/(height(q_UG_site1)/365);
annual_Q2_site1 = sumQ2_site1/(height(q_method2_site1)/365);
annual_Q3_site1 = sumQ3_site1/(height(q_method3_site1)/365);
annual_Q4_site1 = sumQ4_site1/(height(q_method4_site1)/365);

%% unit gradient approach for site 2 (Treated woodland) using sensor MP data
daily_mean_MP_site2_cm = readtimetable("daily_mean_MP_site2_cm.xlsx");
% For site2
theta_r = VGparameter(6:10,1);
theta_s = VGparameter(6:10,2);
alpha = VGparameter(6:10,3);
n = VGparameter(6:10,4);
Ks = VGparameter(6:10,5); % cm/day
m_VG = 1-1./n;
L=0.5;

% Calculate K
K_psi_store = zeros(height(daily_mean_MP_site2_cm),5);

% try K(Se), when L= 0.5, same to the K(h)
for i = 1:5
    alpha_psi = abs(alpha(i).*table2array(daily_mean_MP_site2_cm(:,i))); % to solve the complex double results, should be absolute value!
%     S = (1/(1+alpha_psi.^n(i))).^m_VG(i);
    S = (1+(alpha_psi).^n(i)).^(-m_VG(i));
    K_psi = Ks(i).*S.^L.*(1-(1-S.^(1./m_VG(i))).^m_VG(i)).^2;
    % Store the K values
    K_psi_store(:,i) = K_psi; % cm/day
end

q_UG_site2 = -K_psi_store(:,5); % K at 100 cm (q=-K()) -- cm/day
%% Arithmetic mean K for site 2 (treated woodland) using sensor MP data
K_mean_store = zeros(height(daily_mean_MP_site2_cm),4);
for i = 1:4
    K_mean = mean(K_psi_store(:,i:5),2);
    K_mean_store(:,i) = K_mean;
end

h_gradient_20_100 = (daily_mean_MP_site2_cm.corrected_MP_site2_cm5 - daily_mean_MP_site2_cm.corrected_MP_site2_cm1)/(-80);
h_gradient_40_100 = (daily_mean_MP_site2_cm.corrected_MP_site2_cm5 - daily_mean_MP_site2_cm.corrected_MP_site2_cm2)/(-60);
h_gradient_60_100 = (daily_mean_MP_site2_cm.corrected_MP_site2_cm5 - daily_mean_MP_site2_cm.corrected_MP_site2_cm3)/(-40);
h_gradient_80_100 = (daily_mean_MP_site2_cm.corrected_MP_site2_cm5 - daily_mean_MP_site2_cm.corrected_MP_site2_cm4)/(-20);
h_gradient = [h_gradient_20_100 h_gradient_40_100 h_gradient_60_100 h_gradient_80_100];

q_method2_site2 = -K_mean_store.*(1+h_gradient);
%% Layer-based K with harmonic mean for site 2 (treated woodland) using sensor MP data
layerK_store = zeros(height(daily_mean_MP_site2_cm),4);
for i = 1:4
    layerK = mean(K_psi_store(:,[i,i+1]),2);
    layerK_store(:,i) = layerK;
end

K_eff_store = zeros(height(daily_mean_MP_site2_cm),3);
for i = 1:3
    K_eff = harmmean(layerK_store(:,i:4),2);
    K_eff_store(:,i) = K_eff;
end

q_method3_site2 = -K_eff_store.*(1+h_gradient(:,1:3));%cm/day

% Create a new timetable combining the time and q_UG_site4
q_combined_site2 = timetable(daily_mean_MP_site2_cm.Time, q_UG_site2,q_method2_site2,q_method3_site2);
% Rename the variable for clarity (optional)
q_combined_site2.Properties.VariableNames = {'q_UG','q_method2','q_method3'};

sumQ1_site2 = sum(q_UG_site2);
sumQ2_site2 = sum(q_method2_site2);
sumQ3_site2 = sum(q_method3_site2);
% if the above drainage rates are too high, may need to check the VG
% parameters/or check the drainage rate calculation (+/-)
annual_Q1_site2 = sumQ1_site2/(height(q_UG_site2)/365);
annual_Q2_site2 = sumQ2_site2/(height(q_method2_site2)/365);
annual_Q3_site2 = sumQ3_site2/(height(q_method3_site2)/365);
%% K with harmonic mean for site 2 using sensor MP data
K_eff_store = zeros(height(daily_mean_MP_site2_cm),4);
% Compute harmonic mean for depths 20-100 cm, 40-100 cm, and 60-100 cm
for i = 1:4
    K_eff = harmmean(K_psi_store(:,i:5),2);
    K_eff_store(:,i) = K_eff;
end

q_method4_site2 = -K_eff_store.*(1+h_gradient);%cm/day

q_combined_site2 = q_combined_site2(1:796,:);
% Create a new timetable combining the time and q_UG_site4
q_combined_site2 = timetable(daily_mean_MP_site2_cm.Time, q_UG_site2,q_method2_site2,q_method3_site2,q_method4_site2);
% Rename the variable for clarity (optional)
q_combined_site2.Properties.VariableNames = {'q_UG','q_method2','q_method3','q_method4'};

sumQ1_site2 = sum(q_UG_site2);
sumQ2_site2 = sum(q_method2_site2);
sumQ3_site2 = sum(q_method3_site2);
sumQ4_site2 = sum(q_method4_site2);

annual_Q1_site2 = sumQ1_site2/(height(q_UG_site2)/365);
annual_Q2_site2 = sumQ2_site2/(height(q_method2_site2)/365);
annual_Q3_site2 = sumQ3_site2/(height(q_method3_site2)/365);
annual_Q4_site2 = sumQ4_site2/(height(q_method4_site2)/365);

%% unit gradient approach for site 3 using sensor MP data
daily_mean_MP_site3_cm = readtimetable("daily_mean_MP_site3_cm.xlsx");
% For site3
theta_r = VGparameter(11:15,1);
theta_s = VGparameter(11:15,2);
alpha = VGparameter(11:15,3);
n = VGparameter(11:15,4);
Ks = VGparameter(11:15,5); % cm/day
m_VG = 1-1./n;
L=0.5;

% Calculate K
K_psi_store = zeros(height(daily_mean_MP_site3_cm),5);

% try K(Se), when L= 0.5, same to the K(h)
for i = 1:5
    alpha_psi = abs(alpha(i).*table2array(daily_mean_MP_site3_cm(:,i))); % to solve the complex double results, should be absolute value!
%     S = (1/(1+alpha_psi.^n(i))).^m_VG(i);
    S = (1+(alpha_psi).^n(i)).^(-m_VG(i));
    K_psi = Ks(i).*S.^L.*(1-(1-S.^(1./m_VG(i))).^m_VG(i)).^2;
    % Store the K values
    K_psi_store(:,i) = K_psi; % cm/day
end

q_UG_site3 = -K_psi_store(:,5); % K at 100 cm (q=-K()) -- cm/day

%% Arithmetic mean K for site 3 using sensor MP data
K_mean_store = zeros(height(daily_mean_MP_site3_cm),4);
for i = 1:4
    K_mean = mean(K_psi_store(:,i:5),2);
    K_mean_store(:,i) = K_mean;
end

h_gradient_20_100 = (daily_mean_MP_site3_cm.corrected_MP_site3_cm5 - daily_mean_MP_site3_cm.corrected_MP_site3_cm1)/(-80);
h_gradient_40_100 = (daily_mean_MP_site3_cm.corrected_MP_site3_cm5 - daily_mean_MP_site3_cm.corrected_MP_site3_cm2)/(-60);
h_gradient_60_100 = (daily_mean_MP_site3_cm.corrected_MP_site3_cm5 - daily_mean_MP_site3_cm.corrected_MP_site3_cm3)/(-40);
h_gradient_80_100 = (daily_mean_MP_site3_cm.corrected_MP_site3_cm5 - daily_mean_MP_site3_cm.corrected_MP_site3_cm4)/(-20);
h_gradient = [h_gradient_20_100 h_gradient_40_100 h_gradient_60_100 h_gradient_80_100];

q_method2_site3 = -K_mean_store.*(1+h_gradient);
%% Layer-based K with harmonic mean for site 3 using sensor MP data
layerK_store = zeros(height(daily_mean_MP_site3_cm),4);
for i = 1:4
    layerK = mean(K_psi_store(:,[i,i+1]),2);
    layerK_store(:,i) = layerK;
end

K_eff_store = zeros(height(daily_mean_MP_site3_cm),3);
for i = 1:3
    K_eff = harmmean(layerK_store(:,i:4),2);
    K_eff_store(:,i) = K_eff;
end

q_method3_site3 = -K_eff_store.*(1+h_gradient(:,1:3));%cm/day

% Create a new timetable combining the time and q_UG_site4
q_combined_site3 = timetable(daily_mean_MP_site3_cm.Time, q_UG_site3,q_method2_site3,q_method3_site3);
% Rename the variable for clarity (optional)
q_combined_site3.Properties.VariableNames = {'q_UG','q_method2','q_method3'};

sumQ1_site3 = sum(q_UG_site3);
sumQ2_site3 = sum(q_method2_site3);
sumQ3_site3 = sum(q_method3_site3);

annual_Q1_site3 = sumQ1_site3/(height(q_UG_site3)/365);
annual_Q2_site3 = sumQ2_site3/(height(q_method2_site3)/365);
annual_Q3_site3 = sumQ3_site3/(height(q_method3_site3)/365);
%% K with harmonic mean for site 3 using sensor MP data
K_eff_store = zeros(height(daily_mean_MP_site3_cm),4);
% Compute harmonic mean for depths 20-100 cm, 40-100 cm, and 60-100 cm
for i = 1:4
    K_eff = harmmean(K_psi_store(:,i:5),2);
    K_eff_store(:,i) = K_eff;
end

q_method4_site3 = -K_eff_store.*(1+h_gradient);%cm/day

% Create a new timetable combining the time and q_UG_site4
q_combined_site3 = timetable(daily_mean_MP_site3_cm.Time, q_UG_site3,q_method2_site3,q_method3_site3,q_method4_site3);
% Rename the variable for clarity (optional)
q_combined_site3.Properties.VariableNames = {'q_UG','q_method2','q_method3','q_method4'};

sumQ1_site3 = sum(q_UG_site3);
sumQ2_site3 = sum(q_method2_site3);
sumQ3_site3 = sum(q_method3_site3);
sumQ4_site3 = sum(q_method4_site3);

annual_Q1_site3 = sumQ1_site3/(height(q_UG_site3)/365);
annual_Q2_site3 = sumQ2_site3/(height(q_method2_site3)/365);
annual_Q3_site3 = sumQ3_site3/(height(q_method3_site3)/365);
annual_Q4_site3 = sumQ4_site3/(height(q_method4_site3)/365);

%% unit gradient approach for site 4 using sensor MP data
daily_mean_MP_site4_cm = readtimetable("daily_mean_MP_site4_cm.xlsx");
% For site4
theta_r = VGparameter(16:20,1);
theta_s = VGparameter(16:20,2);
alpha = VGparameter(16:20,3);
n = VGparameter(16:20,4);
Ks = VGparameter(16:20,5); % cm/day
m_VG = 1-1./n;
L=0.5;

% Calculate K
K_psi_store = zeros(height(daily_mean_MP_site4_cm),5);

% try K(Se), when L= 0.5, same to the K(h)
for i = 1:5
    alpha_psi = abs(alpha(i).*table2array(daily_mean_MP_site4_cm(:,i))); % to solve the complex double results, should be absolute value!
%     S = (1/(1+alpha_psi.^n(i))).^m_VG(i);
    S = (1+(alpha_psi).^n(i)).^(-m_VG(i));
    K_psi = Ks(i).*S.^L.*(1-(1-S.^(1./m_VG(i))).^m_VG(i)).^2;
    % Store the K values
    K_psi_store(:,i) = K_psi; % cm/day
end

q_UG_site4 = -K_psi_store(:,5); % K at 100 cm (q=-K()) -- cm/day
%% Arithmetic mean K for site 4 using sensor MP data
K_mean_store = zeros(height(daily_mean_MP_site4_cm),4);
for i = 1:4
    K_mean = mean(K_psi_store(:,i:5),2);
    K_mean_store(:,i) = K_mean;
end

h_gradient_20_100 = (daily_mean_MP_site4_cm.corrected_MP_site4_cm5 - daily_mean_MP_site4_cm.corrected_MP_site4_cm1)/(-80);
h_gradient_40_100 = (daily_mean_MP_site4_cm.corrected_MP_site4_cm5 - daily_mean_MP_site4_cm.corrected_MP_site4_cm2)/(-60);
h_gradient_60_100 = (daily_mean_MP_site4_cm.corrected_MP_site4_cm5 - daily_mean_MP_site4_cm.corrected_MP_site4_cm3)/(-40);
h_gradient_80_100 = (daily_mean_MP_site4_cm.corrected_MP_site4_cm5 - daily_mean_MP_site4_cm.corrected_MP_site4_cm4)/(-20);
h_gradient = [h_gradient_20_100 h_gradient_40_100 h_gradient_60_100 h_gradient_80_100];

q_method2_site4 = -K_mean_store.*(1+h_gradient);
%% Layer-based K with harmonic mean for site 4 using sensor MP data
layerK_store = zeros(height(daily_mean_MP_site4_cm),4);
for i = 1:4
    layerK = mean(K_psi_store(:,[i,i+1]),2);
    layerK_store(:,i) = layerK;
end

K_eff_store = zeros(height(daily_mean_MP_site4_cm),3);
for i = 1:3
    K_eff = harmmean(layerK_store(:,i:4),2);
    K_eff_store(:,i) = K_eff;
end

q_method3_site4 = -K_eff_store.*(1+h_gradient(:,1:3));%cm/day

% Create a new timetable combining the time and q_UG_site4
q_combined_site4 = timetable(daily_mean_MP_site4_cm.Time, q_UG_site4,q_method2_site4,q_method3_site4);
% Rename the variable for clarity (optional)
q_combined_site4.Properties.VariableNames = {'q_UG','q_method2','q_method3'};

sumQ1_site4 = sum(q_UG_site4);
sumQ2_site4 = sum(q_method2_site4);
sumQ3_site4 = sum(q_method3_site4);

annual_Q1_site4 = sumQ1_site4/(height(q_UG_site4)/365);
annual_Q2_site4 = sumQ2_site4/(height(q_method2_site4)/365);
annual_Q3_site4 = sumQ3_site4/(height(q_method3_site4)/365);
%% K with harmonic mean for site 4 using sensor MP data
K_eff_store = zeros(height(daily_mean_MP_site4_cm),4);
% Compute harmonic mean for depths 20-100 cm, 40-100 cm, and 60-100 cm
for i = 1:4
    K_eff = harmmean(K_psi_store(:,i:5),2);
    K_eff_store(:,i) = K_eff;
end

q_method4_site4 = -K_eff_store.*(1+h_gradient);%cm/day

% Create a new timetable combining the time and q_UG_site4
q_combined_site4 = timetable(daily_mean_MP_site4_cm.Time, q_UG_site4,q_method2_site4,q_method3_site4,q_method4_site4);
% Rename the variable for clarity (optional)
q_combined_site4.Properties.VariableNames = {'q_UG','q_method2','q_method3','q_method4'};

sumQ1_site4 = sum(q_UG_site4);
sumQ2_site4 = sum(q_method2_site4);
sumQ3_site4 = sum(q_method3_site4);
sumQ4_site4 = sum(q_method4_site4);

annual_Q1_site4 = sumQ1_site4/(height(q_UG_site4)/365);
annual_Q2_site4 = sumQ2_site4/(height(q_method2_site4)/365);
annual_Q3_site4 = sumQ3_site4/(height(q_method3_site4)/365);
annual_Q4_site4 = sumQ4_site4/(height(q_method4_site4)/365);

%% fill in the missing q for site 1 and site 2 from site 3&4
subsQ_site3 = q_combined_site3(472:end,:);
subsQ_site4 = q_combined_site4(472:end,:);
% Compute element-wise mean
subsQ_avg_site3_site4 = subsQ_site3; % Retain the timetable structure
subsQ_avg_site3_site4{:,:}= (subsQ_site3{:,:} + subsQ_site4{:,:}) / 2;

% subsitute mean of site 3&4 to site 2 (09/04/24-09/30/24)
% Extract the subset from 09/04/24 - 09/30/24
subsQ_avg4site2 = subsQ_avg_site3_site4(21:end, :);
% Add substitute Q for site 2
q_combined_site2 = q_combined_site2(1:796,:);
q_combined_site2_full = vertcat(q_combined_site2,subsQ_avg4site2);

% extract the date (08/15-09/30)
subsQ_site2 = q_combined_site2_full(777:end,:);
subsQ_avg4site1 = subsQ_site3; % Retain the timetable structure
subsQ_avg4site1{:,:} = (subsQ_site3{:,:} + subsQ_site4{:,:} + subsQ_site2{:,:}) / 3;
% add substitute Q for site 1
q_combined_site1_full = vertcat(q_combined_site1,subsQ_avg4site1);

%% plot MP at each depths across 4 sites with rainfall
depths = [20, 40, 60, 80, 100];

% Store the timetables in a cell array
timetables = {daily_mean_MP_site1_cm, daily_mean_MP_site2_cm, daily_mean_MP_site3_cm, daily_mean_MP_site4_cm};

% Site names for display
siteNames = {'Site 1 - Savanna grassland', 'Site 2 - thicketized woodland 1', 'Site 3 - open woodland', 'Site 4 - thicketized woodland 2'};

% Loop through each depth and plot the data from all four timetables
for i = 1:length(depths)
    figure;
    hold on
    
    % Plot matric potential data for all sites
    for j = 1:length(timetables)
        time = timetables{j}.Time;
        data = table2array(timetables{j}(:, i));
        
        % Plot matric potential for the current depth from each site
        plot(time, data, 'DisplayName', [siteNames{j}], 'LineWidth', 1.5);
    end
    
    % Set the left y-axis for matric potential
    xlabel('Time');
    ylabel(['Matric potential (cm) at ' num2str(depths(i)) ' cm depth']);
    set(gca, 'YScale', 'log');
    % Plot rainfall time series on the secondary y-axis
    daily_rain = readtable('daily_rain_CR6.csv');
    yyaxis right;  % Activate the right y-axis
    bar(daily_rain.TIMESTAMP, daily_rain.Rain_cm_Tot, 'c', 'BarWidth', 1.5, 'DisplayName', 'Rainfall (cm)','FaceAlpha', 0.8);
    % Reverse the y-axis for rainfall
    set(gca, 'YDir', 'reverse');
    % Set the color of the right y-axis to black (or any color of your choice)
    ax = gca;
    ax.YAxis(2).Color = 'k';  % Change the secondary (right) y-axis color to black
    ylabel('Rainfall (cm)');  % Label the right y-axis
    
    % Add a combined legend for both matric potential and rainfall
    legend_handle = legend('Location', 'best', 'FontSize', 10);
    legend_handle.ItemTokenSize = [5, 5];  % Adjust line length in the legend
    
    title(['Comparison of Matric Potential and Rainfall at ' num2str(depths(i)) ' cm Depth Across 4 Sites']);
    
    hold off;
end

%% Calculate drainage rates per water year (WY23,WY24)
% Define the start and end dates for Water Years (WY23, WY24)
WY23_start = datetime(2022, 10, 1);  % Start of WY23
WY23_end = datetime(2023, 9, 30);    % End of WY23
WY24_start = datetime(2023, 10, 1);  % Start of WY24
WY24_end = datetime(2024, 9, 30);    % End of WY24

% List of timetables for each site
sites = {'q_combined_site1_full', 'q_combined_site2_full', 'q_combined_site3', 'q_combined_site4'};
site_names = {'Site1', 'Site2', 'Site3', 'Site4'};
% Define the sites to include for WY23 calculation
sites_WY23 = {'q_combined_site1_full', 'q_combined_site2_full'};
site_names_WY23 = {'Site1', 'Site2'};

% Initialize arrays to store the results
sum_WY23_sites12 = zeros(length(sites_WY23), 12);  % 4x3 for each site, 3 methods
sum_WY24_all_sites = zeros(length(sites), 12);
WY23_q_method4_20 = zeros(365,2);
WY24_q_method4_20 = zeros(366,4);

% Loop through each site for WY23 (Site 1 and Site 2 only)
for s = 1:length(sites_WY23)
    site_data = eval(sites_WY23{s});  % Get the timetable for the current site
    
    % Filter the timetable for Water Year 2023 (WY23)
    WY23_data = site_data(site_data.Time >= WY23_start & site_data.Time <= WY23_end, :);

    % Filter the timetable for Water Year 2023 (WY23) -- q_method4_20
    WY23_q_method4 = site_data(site_data.Time >= WY23_start & site_data.Time <= WY23_end, "q_method4");
    WY23_q_method4_20(:,s) =  WY23_q_method4.q_method4(:,1); % 1st column is site 1; 2nd column is site 2
    
    % Calculate the sum of drainage rates for WY23
    sum_WY23 = varfun(@sum, WY23_data, 'InputVariables', {'q_UG', 'q_method2', 'q_method3', 'q_method4'});
    sum_WY23_sites12(s, :) = sum_WY23.Variables;  % Store the result
end
% Convert the results for WY23 into a table with meaningful column names for Site 1 and Site 2 only
WY23_table_sites12 = table(site_names_WY23', sum_WY23_sites12(:,1), sum_WY23_sites12(:,2), sum_WY23_sites12(:,3), sum_WY23_sites12(:,4), ...
                           sum_WY23_sites12(:,5), sum_WY23_sites12(:,6), sum_WY23_sites12(:,7), sum_WY23_sites12(:,8), ...
                           sum_WY23_sites12(:,9), sum_WY23_sites12(:,10), sum_WY23_sites12(:,11), sum_WY23_sites12(:,12), ...
                           'VariableNames', {'Site', 'q_UG', 'q_method2_20', 'q_method2_40', 'q_method2_60', 'q_method2_80', ...
                                             'q_method3_20', 'q_method3_40', 'q_method3_60', ...
                                             'q_method4_20', 'q_method4_40', 'q_method4_60', 'q_method4_80'});

% Loop through each site for WY24
for s = 1:length(sites)
    site_data = eval(sites{s});  % Get the timetable for the current site
      
    % Filter the timetable for Water Year 2024 (WY24)
    WY24_data = site_data(site_data.Time >= WY24_start & site_data.Time <= WY24_end, :);

    % Filter the timetable for Water Year 2024 (WY24) -- q_method4_20
    WY24_q_method4 = site_data(site_data.Time >= WY24_start & site_data.Time <= WY24_end, "q_method4");
    WY24_q_method4_20(:,s) =  WY24_q_method4.q_method4(:,1); % 1st/2nd/3rd/4th column is site 1/2/3/4
    
    % Calculate the sum of drainage rates for WY24
    sum_WY24 = varfun(@sum, WY24_data, 'InputVariables', {'q_UG', 'q_method2', 'q_method3','q_method4'});
    sum_WY24_all_sites(s, :) = sum_WY24.Variables;  % Store the result
end

% Convert the results for WY24 into a table with the same column names
WY24_table = table(site_names', sum_WY24_all_sites(:,1), sum_WY24_all_sites(:,2), sum_WY24_all_sites(:,3), sum_WY24_all_sites(:,4), ...
                   sum_WY24_all_sites(:,5), sum_WY24_all_sites(:,6), sum_WY24_all_sites(:,7), sum_WY24_all_sites(:,8), ...
                   sum_WY24_all_sites(:,9), sum_WY24_all_sites(:,10), sum_WY24_all_sites(:,11), sum_WY24_all_sites(:,12), ...
                   'VariableNames', {'Site', 'q_UG', 'q_method2_20', 'q_method2_40', 'q_method2_60', 'q_method2_80', ...
                                     'q_method3_20', 'q_method3_40', 'q_method3_60',...
                                     'q_method4_20', 'q_method4_40', 'q_method4_60', 'q_method4_80'});

% Add a new column to distinguish between WY23 and WY24
WY23_table_sites12.WaterYear = repmat({'WY23'}, height(WY23_table_sites12), 1);
WY24_table.WaterYear = repmat({'WY24'}, height(WY24_table), 1);

% Concatenate the two tables vertically
final_table = [WY23_table_sites12; WY24_table];

% Reorder the columns to have 'WaterYear' as the first column
final_table = final_table(:, [end, 1, 2:13]);

%% How much rainfall contribute to drainage
% Load rainfall data
daily_rain = readtable('daily_rain_CR6.csv');  % Ensure the table has 'Date' and 'Rainfall' columns

% Filter the rainfall data for each water year
rainfall_WY23 = daily_rain(daily_rain.TIMESTAMP >= WY23_start & daily_rain.TIMESTAMP <= WY23_end, :);
rainfall_WY24 = daily_rain(daily_rain.TIMESTAMP >= WY24_start & daily_rain.TIMESTAMP <= WY24_end, :);

% Sum the rainfall for each water year
total_rainfall_WY23 = sum(rainfall_WY23.Rain_cm_Tot);  % Total rainfall in WY23
total_rainfall_WY24 = sum(rainfall_WY24.Rain_cm_Tot);  % Total rainfall in WY24

% Calculate the contribution of rainfall to drainage
% for WY23 (Site 1 and Site 2)
P2q_WY23_sites12 = abs(sum_WY23_sites12) / total_rainfall_WY23 * 100;  % WY23 contribution ratio
% for WY24
P2q_WY24 = abs(sum_WY24_all_sites) / total_rainfall_WY24 *100;  % WY24 contribution ratio

% Create tables for WY23 and WY24
% Convert the results for WY23 into a table with meaningful column names
P2q_WY23_table_sites12 = table(site_names_WY23', P2q_WY23_sites12(:,1), P2q_WY23_sites12(:,2), P2q_WY23_sites12(:,3), P2q_WY23_sites12(:,4),...
    P2q_WY23_sites12(:,5), P2q_WY23_sites12(:,6),P2q_WY23_sites12(:,7),P2q_WY23_sites12(:,8), ...
                        P2q_WY23_sites12(:,9),P2q_WY23_sites12(:,10), P2q_WY23_sites12(:,11), P2q_WY23_sites12(:,12), ...
                   'VariableNames', {'Site', 'P2q_q_UG', 'P2q_method2_20', 'P2q_method2_40', 'P2q_method2_60', 'P2q_method2_80', ...
                                     'P2q_method3_20', 'P2q_method3_40', 'P2q_method3_60',...
                                     'P2q_method4_20', 'P2q_method4_40', 'P2q_method4_60', 'P2q_method4_80'});

% Convert the results for WY24 into a table with the same column names
P2q_WY24_table = table(site_names', P2q_WY24(:,1), P2q_WY24(:,2), P2q_WY24(:,3), P2q_WY24(:,4), P2q_WY24(:,5), P2q_WY24(:,6), P2q_WY24(:,7), P2q_WY24(:,8), ...
                        P2q_WY24(:,9), P2q_WY24(:,10), P2q_WY24(:,11), P2q_WY24(:,12), ...
                   'VariableNames', {'Site', 'P2q_q_UG', 'P2q_method2_20', 'P2q_method2_40', 'P2q_method2_60', 'P2q_method2_80', ...
                                     'P2q_method3_20', 'P2q_method3_40', 'P2q_method3_60',...
                                     'P2q_method4_20', 'P2q_method4_40', 'P2q_method4_60', 'P2q_method4_80'});

% Add a new column to distinguish between WY23 and WY24
P2q_WY23_table_sites12.WaterYear = repmat({'WY23'}, height(P2q_WY23_table_sites12), 1);
P2q_WY24_table.WaterYear = repmat({'WY24'}, height(P2q_WY24_table), 1);

% Concatenate the two tables vertically
P2q_final_table = [P2q_WY23_table_sites12; P2q_WY24_table];

% Reorder the columns to have 'WaterYear' as the first column
P2q_final_table = P2q_final_table(:, [end, 1, 2:13]);

%% Neutron probe based calculation ########################################
%% VWC from neutron probe on Site A (savanna grassland)
neutronR = readtable("NMM measurements_GusEngeling_matlab.xlsx");

% Convert the cell array to a matrix of doubles
neutronR{:, 19:27} = cellfun(@str2double, neutronR{:, 19:27}, 'UniformOutput', false);
newtable = cell2table(neutronR{:, 19:27});
newtable.Properties.VariableNames = ["Var19","Var20","Var21","Var22","Var23","Var24","Var25","Var26","Var27"];
Orig = neutronR(:,1:18);
neutronR = [Orig,newtable];

%% 
% separate data records for each site
A1_R = neutronR(strcmp(neutronR.("Var4"), 'A1'), :);
A2_R = neutronR(strcmp(neutronR.("Var4"), 'A2'), :);
A3_R = neutronR(strcmp(neutronR.("Var4"), 'A3'), :);
A4_R = neutronR(strcmp(neutronR.("Var4"), 'A4'), :);
B1_R = neutronR(strcmp(neutronR.("Var4"), 'B1'), :);
B2_R = neutronR(strcmp(neutronR.("Var4"), 'B2'), :);
B3_R = neutronR(strcmp(neutronR.("Var4"), 'B3'), :);
B4_R = neutronR(strcmp(neutronR.("Var4"), 'B4'), :);
B5_R = neutronR(strcmp(neutronR.("Var4"), 'B5'), :);
B6_R = neutronR(strcmp(neutronR.("Var4"), 'B6'), :);
C1_R = neutronR(strcmp(neutronR.("Var4"), 'C1'), :);
C2_R = neutronR(strcmp(neutronR.("Var4"), 'C2'), :);
C3_R = neutronR(strcmp(neutronR.("Var4"), 'C3'), :);

% extract date for each site
A1_time = datetime(table2array(A1_R(:,1)),'InputFormat','dd-MMM-yyyy', 'Format','MM/dd/yy');
A2_time = datetime(table2array(A2_R(:,1)),'InputFormat','dd-MMM-yyyy', 'Format','MM/dd/yy');
A3_time = datetime(table2array(A3_R(:,1)),'InputFormat','dd-MMM-yyyy', 'Format','MM/dd/yy');
A4_time = datetime(table2array(A4_R(:,1)),'InputFormat','dd-MMM-yyyy', 'Format','MM/dd/yy');
B1_time = datetime(table2array(B1_R(:,1)),'InputFormat','dd-MMM-yyyy', 'Format','MM/dd/yy');
B2_time = datetime(table2array(B2_R(:,1)),'InputFormat','dd-MMM-yyyy', 'Format','MM/dd/yy');
B3_time = datetime(table2array(B3_R(:,1)),'InputFormat','dd-MMM-yyyy', 'Format','MM/dd/yy');
B4_time = datetime(table2array(B4_R(:,1)),'InputFormat','dd-MMM-yyyy', 'Format','MM/dd/yy');
B5_time = datetime(table2array(B5_R(:,1)),'InputFormat','dd-MMM-yyyy', 'Format','MM/dd/yy');
B6_time = datetime(table2array(B6_R(:,1)),'InputFormat','dd-MMM-yyyy', 'Format','MM/dd/yy');
C1_time = datetime(table2array(C1_R(:,1)),'InputFormat','dd-MMM-yyyy', 'Format','MM/dd/yy');
C2_time = datetime(table2array(C2_R(:,1)),'InputFormat','dd-MMM-yyyy', 'Format','MM/dd/yy');
C3_time = datetime(table2array(C3_R(:,1)),'InputFormat','dd-MMM-yyyy', 'Format','MM/dd/yy');

% convert neutron count ratio to VWC
% coefficients of neutron probe calibration equation
a = 0.2035; b = 0.0184; 
A1_theta = 0.2035*table2array(A1_R(:, 6:16))-0.0184; % averaged from lowest and middle -- currently 
A2_theta = a*table2array(A2_R(:, 6:21))-b;
A3_theta = a*table2array(A3_R(:, 6:15))-b;
A4_theta = a*table2array(A4_R(:, 6:15))-b;
B1_theta = a*table2array(B1_R(:, 6:19))-b;
B2_theta = a*table2array(B2_R(:, 6:19))-b;
B3_theta = a*table2array(B3_R(:, 6:17))-b;
B4_theta = a*table2array(B4_R(:, 6:15))-b;
B5_theta = a*table2array(B5_R(:, 6:20))-b;
B6_theta = a*table2array(B6_R(:, 6:18))-b;
C1_theta = a*table2array(C1_R(:, 6:18))-b;
C2_theta = a*table2array(C2_R(:, 6:18))-b;
C3_theta = a*table2array(C3_R(:, 6:22))-b;    

% 2 values at 20 cm on A1 are negative, substitute it with minimum positive
min_positive = min(A1_theta(A1_theta > 0));
A1_theta(A1_theta < 0) = min_positive;

% 1 value on C1 is >0.2, likely submerged in the water
C1_theta(C1_theta > 0.2) = NaN;

% combine date with converted VWC for each site
A1_theta_date = [A1_R(:,1),array2table(A1_theta)];
A2_theta_date = [A2_R(:,1),array2table(A2_theta)];
A3_theta_date = [A3_R(:,1),array2table(A3_theta)];
A4_theta_date = [A4_R(:,1),array2table(A4_theta)];
B1_theta_date = [B1_R(:,1),array2table(B1_theta)];
B2_theta_date = [B2_R(:,1),array2table(B2_theta)];
B3_theta_date = [B3_R(:,1),array2table(B3_theta)];
B4_theta_date = [B4_R(:,1),array2table(B4_theta)];
B5_theta_date = [B5_R(:,1),array2table(B5_theta)];
B6_theta_date = [B6_R(:,1),array2table(B6_theta)];
C1_theta_date = [C1_R(:,1),array2table(C1_theta)];
C2_theta_date = [C2_R(:,1),array2table(C2_theta)];
C3_theta_date = [C3_R(:,1),array2table(C3_theta)];

%% Average across site VWC profile plot during each season (1-3,4-6,7-9,10-12) -- years are not differentiated!! not consider variation across years
% Define seasons
seasons = {1:3, 4:6, 7:9, 10:12}; % Seasons as Jan-Mar, Apr-Jun, Jul-Sep, Oct-Dec
seasonLabels = {'Jan-Mar', 'Apr-Jun', 'Jul-Sep', 'Oct-Dec'};

% Exclude unwanted dates
excludedDates = datetime({'10-Nov-2022', '08-Oct-2024','22-Nov-2024'}, 'InputFormat', 'dd-MMM-yyyy');
validIdx = ~ismember(A1_time, excludedDates);
validDates = A1_time(validIdx);

% Get unique years in the dataset
uniqueYears = unique(year(validDates));

% Align data with valid dates
A_sites = {alignWithDates(A1_theta_date, validDates), alignWithDates(A2_theta_date, validDates), ...
           alignWithDates(A3_theta_date, validDates), alignWithDates(A4_theta_date, validDates)};

B_sites = {alignWithDates(B1_theta_date, validDates), alignWithDates(B2_theta_date, validDates), ...
           alignWithDates(B3_theta_date, validDates), alignWithDates(B4_theta_date, validDates), ...
           alignWithDates(B5_theta_date, validDates), alignWithDates(B6_theta_date, validDates)};

C_sites = {alignWithDates(C1_theta_date, validDates), alignWithDates(C2_theta_date, validDates), ...
           alignWithDates(C3_theta_date, validDates)};

% Determine the maximum number of depths across all subsites
maxDepths_A = max(cellfun(@(x) size(x, 2), A_sites));
maxDepths_B = max(cellfun(@(x) size(x, 2), B_sites));
maxDepths_C = max(cellfun(@(x) size(x, 2), C_sites));
maxDepths = max([maxDepths_A, maxDepths_B, maxDepths_C]); % Get global max depth

numSeasons = length(seasons); % Adjust for different years

% Initialize seasonal mean matrices
mean_Aseason = nan(numSeasons, maxDepths);
mean_Bseason = nan(numSeasons, maxDepths);
mean_Cseason = nan(numSeasons, maxDepths);

se_Aseason = nan(numSeasons, maxDepths); % standard error
se_Bseason = nan(numSeasons, maxDepths);
se_Cseason = nan(numSeasons, maxDepths);

% Function to pad each site to maxDepths with NaNs
padToMaxDepth = @(x, maxDepths) [x, nan(size(x, 1), maxDepths - size(x, 2))];

% Loop through each year and each season
seasonCounter = 1;
seasonLabels_new = cell(1, 4);


    for s = 1:length(seasons)
        seasonIdx = ismember(month(validDates), seasons{s}); % Filter for this season in this year

        % Compute seasonal mean for each subsite separately and pad to max depth
        Aseason_means = cellfun(@(x) padToMaxDepth(nanmean(x(seasonIdx, :), 1), maxDepths), A_sites, 'UniformOutput', false); % 1 means computing the mean of the rows (dates)
        Bseason_means = cellfun(@(x) padToMaxDepth(nanmean(x(seasonIdx, :), 1), maxDepths), B_sites, 'UniformOutput', false);
        Cseason_means = cellfun(@(x) padToMaxDepth(nanmean(x(seasonIdx, :), 1), maxDepths), C_sites, 'UniformOutput', false);

        % Convert cell arrays to matrices (now same depth size)
        Aseason_matrix = cell2mat(Aseason_means');
        Bseason_matrix = cell2mat(Bseason_means');
        Cseason_matrix = cell2mat(Cseason_means');

        % Compute final mean across subsites
        mean_Aseason(seasonCounter, :) = nanmean(Aseason_matrix, 1); % for each season, calculate the nanmean of each subsite!!
        mean_Bseason(seasonCounter, :) = nanmean(Bseason_matrix, 1);
        mean_Cseason(seasonCounter, :) = nanmean(Cseason_matrix, 1);

        % Compute standard error -- 0: default normalization (N-1, sample standard deviation)
        se_Aseason(seasonCounter, :) = nanstd(Aseason_matrix, 0, 1) ./ sqrt(sum(~isnan(Aseason_matrix), 1)); % across sites for each depth
        se_Bseason(seasonCounter, :) = nanstd(Bseason_matrix, 0, 1) ./ sqrt(sum(~isnan(Bseason_matrix), 1));
        se_Cseason(seasonCounter, :) = nanstd(Cseason_matrix, 0, 1) ./ sqrt(sum(~isnan(Cseason_matrix), 1));

        % Update season label to include year (e.g., 'Jan-Mar 2023', 'Jan-Mar 2024')
        seasonLabels_new{seasonCounter} = seasonLabels{s};
        
        % Move to next season slot
        seasonCounter = seasonCounter + 1;
    end


% Plot seasonal VWC profiles
depth = 20:20:(20 * maxDepths); % Dynamic depth range

figure;  % Create a new figure for each season
t = tiledlayout(1,4, 'TileSpacing', 'compact', 'Padding', 'compact'); % 4 rows, 1 column
% Define the new subplot order
subplotOrder = [4, 1, 2, 3];
for i = 1:4
    s = subplotOrder(i);  % Get the correct season index
%     subplot(2,2,i);  % Create a 2x2 subplot grid
    nexttile; % Instead of subplot()
   
    % Plot seasonal means with SE as error bars
    h1 = errorbar(mean_Aseason(s, :), depth, se_Aseason(s, :), 'horizontal', '-', 'DisplayName', 'Savanna');
    hold on;
    h2 = errorbar(mean_Bseason(s, :), depth, se_Bseason(s, :), 'horizontal', '-', 'DisplayName', 'Woodland Mosaic');
    h3 = errorbar(mean_Cseason(s, :), depth, se_Cseason(s, :), 'horizontal', '-', 'DisplayName', 'Thicketized Woodland');

    % Store legend handles from the first subplot only
    if i == 1
        legendHandles = [h1, h2, h3];
        ylabel('Depth (cm)');
    end

%     xlabel('θ_v (cm³/cm³)');
    xlabel('\theta_v (cm^3 cm^{-3})');
    
%     title(['VWC Profile - ', seasonLabels_new{s}]);  % Season label now includes year (e.g., 'Jan-Mar 2023')
    title([seasonLabels_new{s}]);  % Season label now includes year (e.g., 'Jan-Mar 2023')
%     legend;
    grid on;
    set(gca, 'XAxisLocation', 'top', 'YDir', 'reverse', 'FontSize', 10);
    axis([0.01 0.16 0 350]);
    xticks(0.01:0.03:0.16);

end

% Create a global legend at the bottom using handles from the first subplot
lg = legend(legendHandles, {'Savanna', 'Woodland Mosaic', 'Thicketized Woodland'}, ...
            'Orientation', 'horizontal', 'Location', 'southoutside');
lg.Layout.Tile = 'south'; % Place legend at the bottom


%% 1. Compute Mean and SE for Each Season & Treatment
% Initialize empty arrays
all_VWC = [];
all_depths = [];
all_treatments = [];
all_seasons = [];

% Loop through seasons and depths
for s = 1:numSeasons
    for d = 1:maxDepths-2
        % Collect VWC values for each treatment
        vwc_savanna = mean_Aseason(s, d);
        vwc_woodland = mean_Bseason(s, d);
        vwc_thicket = mean_Cseason(s, d);

        % Append to dataset (only if not NaN)
        if ~isnan(vwc_savanna)
            all_VWC = [all_VWC; vwc_savanna];
            all_depths = [all_depths; depth(d)];
            all_treatments = [all_treatments; "Savanna"];
            all_seasons = [all_seasons; seasonLabels_new{s}];
        end
        if ~isnan(vwc_woodland)
            all_VWC = [all_VWC; vwc_woodland];
            all_depths = [all_depths; depth(d)];
            all_treatments = [all_treatments; "Woodland Mosaic"];
            all_seasons = [all_seasons; seasonLabels_new{s}];
        end
        if ~isnan(vwc_thicket)
            all_VWC = [all_VWC; vwc_thicket];
            all_depths = [all_depths; depth(d)];
            all_treatments = [all_treatments; "Thicketized Woodland"];
            all_seasons = [all_seasons; seasonLabels_new{s}];
        end
    end
end

% Convert treatment and season to categorical variables
all_treatments = categorical(all_treatments);
all_seasons = categorical(cellstr(all_seasons));
% Compute Mean and Standard Error for Each Treatment in Each Season
% treatments = unique(all_treatments); % Get unique treatment categories
treatments = categorical({'Savanna','Woodland Mosaic','Thicketized Woodland'});
% seasons = unique(all_seasons); % Get unique seasons
seasons = categorical({'Oct-Dec', 'Jan-Mar', 'Apr-Jun', 'Jul-Sep'});

mean_se_table = table(); % Initialize table to store results

for i = 1:length(treatments)
    for j = 1:length(seasons)
        % Extract VWC values for the specific treatment and season
        vwc_values = all_VWC(all_treatments == treatments(i) & all_seasons == seasons(j));
        
        % Compute mean and standard error (SE)
        mean_vwc = mean(vwc_values, 'omitnan'); % Mean
        se_vwc = std(vwc_values, 'omitnan') / sqrt(sum(~isnan(vwc_values))); % Standard Error (SE)
        
        % Append results to the table
        mean_se_table = [mean_se_table; table(treatments(i), seasons(j), mean_vwc, se_vwc)];
    end
end

% Rename table columns
mean_se_table.Properties.VariableNames = {'Treatment', 'Season', 'Mean_VWC', 'SE_VWC'};

%% Run 3-way ANOVA
[p, tbl, stats] = anovan(all_VWC, {all_treatments, all_depths, all_seasons}, ...
    'model', 'interaction', ...
    'varnames', {'Treatment', 'Depth', 'Season'});

% Post-hoc comparison: Treatment × Season
figure;
[c_ts, m_ts, ~, groupnames] = multcompare(stats, 'Dimension', [1, 3]);

% Convert multcompare output to table
comparison_table = array2table(c_ts, ...
    'VariableNames', {'Group1', 'Group2', 'Mean_Diff', 'Lower_CI', 'Upper_CI', 'p_Value'});
comparison_table.Group1_Label = groupnames(c_ts(:,1));
comparison_table.Group2_Label = groupnames(c_ts(:,2));

% Extract Treatment and Season from group labels
label1 = string(comparison_table.Group1_Label);
label2 = string(comparison_table.Group2_Label);
n = height(comparison_table);
T1 = strings(n,1); S1 = strings(n,1);
T2 = strings(n,1); S2 = strings(n,1);

for i = 1:n
    t1 = regexp(label1(i), 'Treatment=(.*?),', 'tokens', 'once');
    s1 = regexp(label1(i), 'Season=(.*)$',    'tokens', 'once');
    t2 = regexp(label2(i), 'Treatment=(.*?),', 'tokens', 'once');
    s2 = regexp(label2(i), 'Season=(.*)$',    'tokens', 'once');

    if ~isempty(t1), T1(i) = t1{1}; end
    if ~isempty(s1), S1(i) = s1{1}; end
    if ~isempty(t2), T2(i) = t2{1}; end
    if ~isempty(s2), S2(i) = s2{1}; end
end

% Add parsed values to table
comparison_table.Treatment1 = T1;
comparison_table.Season1    = S1;
comparison_table.Treatment2 = T2;
comparison_table.Season2    = S2;

% Logical filters
is_significant   = comparison_table.p_Value < 0.05;

same_season      = comparison_table.Season1 == comparison_table.Season2;
diff_treatment   = comparison_table.Treatment1 ~= comparison_table.Treatment2;

same_treatment   = comparison_table.Treatment1 == comparison_table.Treatment2;
diff_season      = comparison_table.Season1 ~= comparison_table.Season2;

% --- 1. Differences between treatments in each season
significant_treatment_diff_by_season = comparison_table(...
    is_significant & same_season & diff_treatment, :);

% --- 2. Differences between seasons within each treatment
significant_season_diff_by_treatment = comparison_table(...
    is_significant & same_treatment & diff_season, :);

%% check ANOVA for depth
[p, tbl, stats] = anovan(all_VWC, {all_treatments, all_depths, all_seasons}, ...
    'model', 'interaction', ...
    'varnames', {'Treatment', 'Depth', 'Season'});
% check p value for each treatment at each depth
figure;
[c_td, m_td, ~, groupnames_td] = multcompare(stats, 'Dimension', [1, 2]);
comparison_td = array2table(c_td, ...
    'VariableNames', {'Group1', 'Group2', 'Mean_Diff', 'Lower_CI', 'Upper_CI', 'p_Value'});
comparison_td.Group1_Label = groupnames_td(c_td(:,1));
comparison_td.Group2_Label = groupnames_td(c_td(:,2));

depths = 20:20:300;   % depths used
group_id = (1:length(depths)*3)';  % 45 total groups
depth_per_group = repelem(depths, 3)';  % repeat each depth 3 times

group_depth_map = table(group_id, depth_per_group, 'VariableNames', {'Group', 'Depth'});

comparison_td = outerjoin(comparison_td, group_depth_map, ...
    'LeftKeys', 'Group1', 'RightKeys', 'Group', 'RightVariables', 'Depth', ...
    'Type', 'left', 'MergeKeys', true);
comparison_td.Properties.VariableNames{'Depth'} = 'Depth1';

comparison_td = outerjoin(comparison_td, group_depth_map, ...
    'LeftKeys', 'Group2', 'RightKeys', 'Group', 'RightVariables', 'Depth', ...
    'Type', 'left', 'MergeKeys', true);
comparison_td.Properties.VariableNames{'Depth'} = 'Depth2';

same_depth = comparison_td.Depth1 == comparison_td.Depth2;
significant = comparison_td.p_Value < 0.05;
different_groups = comparison_td.Group1_Group ~= comparison_td.Group2_Group;

significant_treatment_diff_by_depth = comparison_td(same_depth & different_groups & significant, :);
