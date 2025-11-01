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

% if the above drainage rates are too high, may need to check the VG
% parameters/or check the drainage rate calculation (+/-)

%% unit gradient approach for site 2 (woodland mosaic_under-canopy1) using sensor MP data
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
%% Arithmetic mean K for site 2 (woodland mosaic_under-canopy1) using sensor MP data
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
%% Layer-based K with harmonic mean for site 2 (woodland mosaic_under-canopy1) using sensor MP data
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

%% K with harmonic mean for site 2 (woodland mosaic_under-canopy1) using sensor MP data
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

%% unit gradient approach for site 3 (woodland mosaic_inter-canopy) using sensor MP data
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

%% Arithmetic mean K for site 3 (woodland mosaic_inter-canopy) using sensor MP data
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
%% Layer-based K with harmonic mean for site 3 (woodland mosaic_inter-canopy) using sensor MP data
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

%% K with harmonic mean for site 3 (woodland mosaic_inter-canopy) using sensor MP data
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

%% unit gradient approach for site 4 (woodland mosaic_under-canopy2) using sensor MP data
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
%% Arithmetic mean K for site 4 (woodland mosaic_under-canopy2) using sensor MP data
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
%% Layer-based K with harmonic mean for site 4 (woodland mosaic_under-canopy2) using sensor MP data
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

%% K with harmonic mean for site 4 (woodland mosaic_under-canopy2) using sensor MP data
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

%% plot drainage rate calculated from method4_20: figure 4
figure
plot(q_combined_site1_full.Time(129:end),q_combined_site1_full.q_method4(129:end,1),'LineWidth',1);hold on;
plot(q_combined_site2_full.Time(93:end),q_combined_site2_full.q_method4(93:end,1),'LineWidth',1);hold on;
plot(q_combined_site3.Time(153:end),q_method4_site3(153:end,1),'LineWidth',1);hold on;
plot(q_combined_site4.Time(153:end),q_method4_site4(153:end,1),'LineWidth',1);
legend('Savanna', 'Woodland mosaic 1: under canopy', 'Woodland mosaic 2: intercanopy', 'Woodland mosaic 3: under canopy','location','best');

% Generate seasonal tick marks (every three months)
tick_dates = datetime(2022,10,1):calmonths(6):datetime(2024,10,01); % Adjust x_data accordingly

% Format the x-axis
xticks(tick_dates);
xticklabels(datestr(tick_dates, 'mmm yyyy')); % Format as "Mon YYYY"

% Label the axis
xlabel('Date');
ylabel('Drainage rate (cm/day)');

% Set x-axis limits to start from 05/22/24
start_date = datetime(2024, 5, 22);
end_date = max([q_combined_site1_full.Time(end), q_combined_site2_full.Time(end), ...
                q_combined_site3.Time(end), q_combined_site4.Time(end)]); % Auto detect last date

xlim([start_date, end_date]); % Set the x-axis limit

% Generate seasonal tick marks from 05/22/24 onward
tick_dates = start_date:calmonths(1):end_date; 

% Format the x-axis ticks
xticks(tick_dates);
xticklabels(datestr(tick_dates, 'dd mmm yyyy')); % Format as "Mon YYYY"
% xtickangle(45); % Rotate for better readability

% Label the axes
xlabel('Date');
ylabel('Drainage rate (cm day^{-1})');

grid on; % Optional: Add grid for better readability

% Overlay Rainfall Data (Secondary Y-axis)
yyaxis right;
bar(daily_rain.TIMESTAMP, daily_rain.Rain_cm_Tot, 'FaceColor', 'c', 'EdgeColor', 'none', 'BarWidth', 0.8); % Blue bars

ylabel('Rainfall (cm)', 'Color', 'k'); % Right y-axis label
yticks(yticks); % Keep existing yticks
ax = gca; % Get current axis
ax.YColor = 'k'; % Change y-axis ticks and labels to black
ylim([0 max(daily_rain.Rain_cm_Tot) * 0.8]); % Adjust y-axis for visibility


% Add legend for rainfall separately
legend({'Savanna', 'Woodland mosaic 1: under canopy', 'Woodland mosaic 2: intercanopy', ...
        'Woodland mosaic 3: under canopy', 'Rainfall'}, 'Location', 'best');

hold off;

%% drainage with rainfall events: table 4
start_rainfall_period1 = datetime(2024, 5, 22);
end_rainfall_period1 = datetime(2024, 6, 11);
start_rainfall_period2 = datetime(2024, 7, 7);
end_rainfall_period2 = datetime(2024, 7, 24);

start_q_period1 = datetime(2024, 5, 22);
end_q_period1 = datetime(2024, 6, 15);
start_q_period2 = datetime(2024, 7, 19);
end_q_period2 = datetime(2024, 8, 4);

% Extract rainfall data within the specified periods
rain_period1 = daily_rain.Rain_cm_Tot(daily_rain.TIMESTAMP >= start_rainfall_period1 & daily_rain.TIMESTAMP <= end_rainfall_period1);
rain_period2 = daily_rain.Rain_cm_Tot(daily_rain.TIMESTAMP >= start_rainfall_period2 & daily_rain.TIMESTAMP <= end_rainfall_period2);

% Compute total rainfall for each period
total_rainfall_period1 = sum(rain_period1, 'omitnan');
total_rainfall_period2 = sum(rain_period2, 'omitnan');

% Extract and compute total drainage for each site within period 1
total_drainage_site1_period1 = sum(q_combined_site1_full.q_method4(q_combined_site1_full.Time >= start_q_period1 & q_combined_site1_full.Time <= end_q_period1), 'omitnan');
total_drainage_site2_period1 = sum(q_combined_site2_full.q_method4(q_combined_site2_full.Time >= start_q_period1 & q_combined_site2_full.Time <= datetime(2024, 6, 13)), 'omitnan');
total_drainage_site3_period1 = sum(q_method4_site3(q_combined_site3.Time >= start_q_period1 & q_combined_site3.Time <= end_q_period1), 'omitnan');
total_drainage_site4_period1 = sum(q_method4_site4(q_combined_site4.Time >= start_q_period1 & q_combined_site4.Time <= datetime(2024, 6, 17)), 'omitnan');

% Extract and compute total drainage for each site within period 2
total_drainage_site1_period2 = sum(q_combined_site1_full.q_method4(q_combined_site1_full.Time >= start_q_period2 & q_combined_site1_full.Time <= datetime(2024, 8, 5)), 'omitnan');
total_drainage_site2_period2 = sum(q_combined_site2_full.q_method4(q_combined_site2_full.Time >= start_q_period2 & q_combined_site2_full.Time <= end_q_period2), 'omitnan');
total_drainage_site3_period2 = sum(q_method4_site3(q_combined_site3.Time >= start_q_period2 & q_combined_site3.Time <= end_q_period2), 'omitnan');
total_drainage_site4_period2 = sum(q_method4_site4(q_combined_site4.Time >= start_q_period2 & q_combined_site4.Time <= end_q_period2), 'omitnan');

% Display the results
fprintf('Total Rainfall for Period 1 (May 22 - June 11, 2024): %.2f cm\n', total_rainfall_period1);
fprintf('Total Rainfall for Period 2 (July 10 - July 25, 2024): %.2f cm\n\n', total_rainfall_period2);

fprintf('Total Drainage for Period 1 (May 22 - June 11, 2024):\n');
fprintf('   Savanna: %.2f cm\n', total_drainage_site1_period1);
fprintf('   Woodland Mosaic 1 (Under Canopy): %.2f cm\n', total_drainage_site2_period1);
fprintf('   Woodland Mosaic 2 (Intercanopy): %.2f cm\n', total_drainage_site3_period1);
fprintf('   Woodland Mosaic 3 (Under Canopy): %.2f cm\n\n', total_drainage_site4_period1);

fprintf('Total Drainage for Period 2 (July 10 - July 25, 2024):\n');
fprintf('   Savanna: %.2f cm\n', total_drainage_site1_period2);
fprintf('   Woodland Mosaic 1 (Under Canopy): %.2f cm\n', total_drainage_site2_period2);
fprintf('   Woodland Mosaic 2 (Intercanopy): %.2f cm\n', total_drainage_site3_period2);
fprintf('   Woodland Mosaic 3 (Under Canopy): %.2f cm\n', total_drainage_site4_period2);

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

%% Run 3-way ANOVA + Post-hoc comparison: Treatment × Season
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

%% Run 3-way ANOVA + Post-hoc comparison: Treatment × Depth
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

%% deep wells (GW level + deep well Neutron probe) + ET +P
% Load deep wells data
neutronR_DW = readtable('Deep_well.xlsx');
well = readtable('Processed_Water_Level.xlsx'); % Adjust filename

% Extract groundwater data
well6 = well(:, ["Date", "Water_Level_Well6"]); % BCN
well7 = well(:, ["Date", "Water_Level_Well7"]); % BCS
well8 = well(:, ["Date", "Water_Level_Well8"]); % AW
well9 = well(:, ["Date", "Water_Level_Well9"]); % AE

% Convert groundwater level to depth (WTD: Water Table Depth)
well6.WTD = 512.064 - well6.Water_Level_Well6*100;
well7.WTD = 512.064 - well7.Water_Level_Well7*100;
well8.WTD = 728.472 - well8.Water_Level_Well8*100;
well9.WTD = 646.176 - well9.Water_Level_Well9*100;

% Convert groundwater dates to datetime
GW_time = datetime(well6.Date(1:132), 'InputFormat', 'yyyy-MM-dd');

% Extract soil moisture data for each well
AW_R = neutronR_DW(strcmp(neutronR_DW.("Var2"), 'AW'), :);
AE_R = neutronR_DW(strcmp(neutronR_DW.("Var2"), 'AE'), :);
BCN_R = neutronR_DW(strcmp(neutronR_DW.("Var2"), 'BCN'), :);
BCS_R = neutronR_DW(strcmp(neutronR_DW.("Var2"), 'BCS'), :);

% Convert measurement dates
AW_time = datetime(AW_R.Var1, 'InputFormat', 'dd-MMM-yyyy');
AE_time = datetime(AE_R.Var1, 'InputFormat', 'dd-MMM-yyyy');
BCN_time = datetime(BCN_R.Var1, 'InputFormat', 'dd-MMM-yyyy');
BCS_time = datetime(BCS_R.Var1, 'InputFormat', 'dd-MMM-yyyy');

% Convert neutron count ratio to VWC
a = 0.2035; b = 0.0184;
AW_theta = a * table2array(AW_R(:, 3:end)) - b;
AE_theta = a * table2array(AE_R(:, 3:end)) - b;
BCN_theta = a * table2array(BCN_R(:, 3:end)) - b;
BCS_theta = a * table2array(BCS_R(:, 3:end)) - b;

% Define depth levels (assuming 20 cm increments)
depth_DW = (1:size(AW_theta,2)) * 20;

% Find all unique soil moisture measurement dates
all_dates = unique(sort([AW_time; AE_time; BCN_time; BCS_time]));
all_dates_num = datenum(all_dates); % Convert to numeric format

% Assign fixed-width bars for each soil moisture measurement
AW_aligned = NaN(length(all_dates), size(AW_theta,2));
AE_aligned = NaN(length(all_dates), size(AE_theta,2));
BCN_aligned = NaN(length(all_dates), size(BCN_theta,2));
BCS_aligned = NaN(length(all_dates), size(BCS_theta,2));

[~, AW_idx] = ismember(AW_time, all_dates);
[~, AE_idx] = ismember(AE_time, all_dates);
[~, BCN_idx] = ismember(BCN_time, all_dates);
[~, BCS_idx] = ismember(BCS_time, all_dates);

AW_aligned(AW_idx, :) = AW_theta;
AE_aligned(AE_idx, :) = AE_theta;
BCN_aligned(BCN_idx, :) = BCN_theta;
BCS_aligned(BCS_idx, :) = BCS_theta;

% Assign wells and groundwater data
wells = {'AW', 'AE', 'BCN', 'BCS'};
aligned_data = {AW_aligned, AE_aligned, BCN_aligned, BCS_aligned};
GW_series = {well8.WTD(1:132), well9.WTD(1:132), well6.WTD(1:132), well7.WTD(1:132)};

% ET + P 2024
ET_P_grass2024 = readtable("savanna_ET_P_2024.csv"); % unit: mm
ET_P_wood2024 = readtable("woodland_ET_P_2024.csv");

ET_P_grass2024.ET_cm = ET_P_grass2024{:,2}/10;
ET_P_grass2024.P_cm = ET_P_grass2024{:,3}/10;

ET_P_wood2024.ET_cm = ET_P_wood2024{:,2}/10;
ET_P_wood2024.P_cm = ET_P_wood2024{:,3}/10;

ET_P_grass2024 = ET_P_grass2024(67:274, :);
ET_P_wood2024 = ET_P_wood2024(67:274, :);

%% only run the period with all the variables
start_date = datetime('2024-05-22', 'InputFormat', 'yyyy-MM-dd');

% Filter Groundwater Level Data
GW_filter_idx = GW_time >= start_date;
GW_time = GW_time(GW_filter_idx);

% Apply the same filter to groundwater series
GW_series = {well8.WTD(GW_filter_idx), well9.WTD(GW_filter_idx), ...
             well6.WTD(GW_filter_idx), well7.WTD(GW_filter_idx)};

% Filter Soil Moisture Data
all_dates_filter_idx = all_dates >= start_date;
all_dates = all_dates(all_dates_filter_idx);
all_dates_num = datenum(all_dates);

% Apply the filter to soil moisture matrices
AW_aligned = AW_aligned(all_dates_filter_idx, :);
AE_aligned = AE_aligned(all_dates_filter_idx, :);
BCN_aligned = BCN_aligned(all_dates_filter_idx, :);
BCS_aligned = BCS_aligned(all_dates_filter_idx, :);
aligned_data = {AW_aligned, AE_aligned, BCN_aligned, BCS_aligned};

% Filter ET and Precipitation Data
ET_grass_filter_idx = ET_P_grass2024.DateTime >= start_date;
ET_wood_filter_idx = ET_P_wood2024.DateTime >= start_date;

ET_P_grass2024 = ET_P_grass2024(ET_grass_filter_idx, :);
ET_P_wood2024 = ET_P_wood2024(ET_wood_filter_idx, :);

daily_rain = readtable('daily_rain_CR6.csv');
daily_rain_filter_idx = daily_rain.TIMESTAMP >= start_date;
daily_rain_filter = daily_rain(daily_rain_filter_idx,[1,3]); % P: cm
%% --revised figure (figure 5)
figure;
t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:4
%     subplot(2,2,i);
    ax = nexttile;
    hold on;

    for j = 1:length(all_dates)
        if ~all(isnan(aligned_data{i}(j,:)))
            x_center = all_dates_num(j);  % Use center for text label
            x_left = x_center - 5;
            x_right = x_center + 5;
            x_range = [x_left, x_right, x_right, x_left];

            for k = 1:size(depth_DW, 2)-1
                y_top = depth_DW(k);
                y_bottom = depth_DW(k+1);
                y_range = [y_top, y_top, y_bottom, y_bottom];

                color_value = aligned_data{i}(j, k);
                if ~isnan(color_value)
                    patch(x_range, y_range, color_value, 'EdgeColor', 'none');
                end
            end

            % Add date text above each column
            text(x_center, depth_DW(1) - 20, datestr(all_dates(j), 'mm/dd'), ...
                 'Rotation', 0, 'HorizontalAlignment', 'center', 'fontsize',12);

        end
    end

% Define date range
start_date = datetime(2024, 5, 22);
end_date = datetime(2024, 9, 30);

% Convert to datenum
start_num = datenum(start_date);
end_num = datenum(end_date);

% Generate monthly ticks
tick_dates = start_date:calmonths(1):end_date;
tick_nums = datenum(tick_dates);

% Apply to x-axis
% xlim([start_num, end_num]);                        % Corrected: xlim requires numeric input
xticks(tick_nums);
xticklabels(datestr(tick_dates, 'dd mmm yyyy'));
xtickangle(0);
xlabel('Date');

    yyaxis left;
    gw_handle = plot(datenum(GW_time), GW_series{i}, 'k-', 'LineWidth', 2);
    ylabel('Depth (cm)');
    set(gca, 'YDir', 'reverse');
    axis tight;
    if i == 1 || i == 2
        ylim([20 700]);
        yticks(0:100:700);
    end

    if i == 1 || i == 2
        ET_data = ET_P_grass2024;
    else
        ET_data = ET_P_wood2024;
    end

    yyaxis right;
    
    et_handle = plot(datenum(ET_data.DateTime), ET_data.EnsembleET, 'r-', 'LineWidth', 1.5);
    rain_handle = bar(datenum(daily_rain_filter.TIMESTAMP), daily_rain_filter.Rain_cm_Tot, 'c', 'BarWidth', 0.8, 'EdgeColor', 'none');
    rain_handle.EdgeColor = 'none';
    ylim([0 7]);
    ylabel('ET (mm) & P (cm)'); 

    wells = {'Deep well #1 (Savanna)', 'Deep well #2 (Savanna)', ...
             'Deep well #3 (Woodland mosaic)', 'Deep well #4 (Woodland mosaic)'};
    ax = gca;
    ax.YAxis(2).TickLabelInterpreter = 'tex';  % (optional, for formatting)
    ax.YAxis(2).Color = 'k';                   % sets tick label color to black
    ax.YAxis(2).Label.Color = 'k';  % Set right y-axis label color to black

    title(ax, wells{i}, 'Units', 'normalized', 'HorizontalAlignment', 'left', ...
      'Position', [0, 1.05, 0], 'FontWeight', 'bold','fontsize',14);
    hold off;
end

% After the for-loop ends
lgd = legend([gw_handle, et_handle, rain_handle], ...
             {'Groundwater Table Depth (cm)', 'ET (mm)', 'Rainfall (cm)'}, ...
             'Orientation', 'horizontal', 'Location', 'southoutside','FontSize', 12);

% Create a shared colorbar to the far right
colormap(flipud(parula));
cb = colorbar; % Attach to tiledlayout
cb.Label.String = '\theta_v (cm^3 cm^{-3})';
cb.Label.FontSize = 12;
cb.Label.Fontweight = 'bold';
cb.Label.Interpreter = 'tex';  % To ensure proper rendering of θ_v
caxis([0 0.3]);               % Set color range if needed
 
%% soil water storage of AE and BCS
% Deep well #2 (AE)
NP_AE = aligned_data{1,2};                 % [nTime x nLayers]
AE_top = sum(NP_AE(4:7,1:15), 2) * 20/15;     
AE_mid = sum(NP_AE(4:7,16:23),2)*20/8;
AE_bot = sum(NP_AE(4:7,24:26), 2) * 20/3;  

% AE_top_avg = mean(AE_top);
% AE_mid_avg = mean(AE_mid);
% AE_bot_avg = mean(AE_bot);

ratio_top_mid = mean(AE_mid)/mean(AE_top) % 1.80
ratio_top_bot = mean(AE_bot)/mean(AE_top) % 3.12

% Deep well #4 (BCS)
NP_BCS = aligned_data{1,4};                
BCS_top = sum(NP_BCS([4,5,7],1:16), 2) * 20/16;   
BCS_mid = sum(NP_BCS([4,5,7],17:19), 2) * 20/3;   
BCS_bot = sum(NP_BCS([4,5,7],20:22), 2) * 20/3; 

% BCS_top_avg = mean(BCS_top);
% BCS_mid_avg = mean(BCS_mid);
% BCS_bot_avg = mean(BCS_bot);

ratio_top_mid = mean(BCS_mid)/mean(BCS_top) % 2.74
ratio_top_bot = mean(BCS_bot)/mean(BCS_top) % 4.71

%% convert MP into VWC with upper/lower bounds

% Inputs
start_date = datetime(2024,5,22);
end_date   = datetime(2024,8,14);
tr = timerange(start_date,end_date,'closed');

MP1 = daily_mean_MP_site1_cm(tr,:);   % Savanna (site 1), 5 depths
MP2 = daily_mean_MP_site2_cm(tr,:);   % Woodland (site 2), 5 depths
MP2 = MP2(1:85,:); % manually adjust the interested time period

% MP units cm (suction). If negative, abs() is used below.

% Parameters per depth (1..5)
% --- Savanna (site 1)
sav.mean.thr = [0.035 0.035 0.037 0.036 0.046];
sav.mean.ths = [0.412 0.396 0.394 0.394 0.403];
sav.mean.alp = [0.045 0.045 0.044 0.045 0.042];
sav.mean.n   = [2.056 2.198 2.326 2.224 3.403];

sav.lo.thr = [0.013 0.015 0.017 0.016 0.022]; sav.hi.thr = [0.056 0.055 0.056 0.055 0.071];
sav.lo.ths = [0.370 0.361 0.361 0.361 0.366]; sav.hi.ths = [0.454 0.431 0.426 0.428 0.439];
sav.lo.alp = [0.026 0.028 0.029 0.029 0.026]; sav.hi.alp = [0.078 0.071 0.068 0.069 0.067];
sav.lo.n   = [1.587 1.737 1.857 1.768 2.595]; sav.hi.n   = [2.663 2.782 2.914 2.797 4.464];

% --- Woodland (site 2)
wood.mean.thr = [0.034 0.044 0.038 0.038 0.032];
wood.mean.ths = [0.387 0.396 0.409 0.391 0.389];
wood.mean.alp = [0.044 0.043 0.044 0.044 0.045];
wood.mean.n   = [2.002 3.284 2.372 2.480 1.900];

wood.lo.thr = [0.016 0.023 0.016 0.018 0.015]; wood.hi.thr = [0.052 0.065 0.061 0.058 0.050];
wood.lo.ths = [0.356 0.363 0.369 0.360 0.356]; wood.hi.ths = [0.418 0.429 0.449 0.423 0.422];
wood.lo.alp = [0.029 0.028 0.027 0.029 0.029]; wood.hi.alp = [0.067 0.066 0.074 0.067 0.070];
wood.lo.n   = [1.619 2.589 1.827 1.993 1.512]; wood.hi.n   = [2.476 4.166 3.081 3.086 2.388];

% Convert 
[VWC1_mean, VWC1_lo95, VWC1_hi95] = vg_bounds_from_CI(MP1, sav);
[VWC2_mean, VWC2_lo95, VWC2_hi95] = vg_bounds_from_CI(MP2, wood);

%% 0-1m VWC -> 0-1m SWS -> average SWS every 20 cm
layers = 1:5;        % depths 20,40,60,80,100 cm
fac    = 20;         % θ * 20 = cm water per layer

% ---------- Savanna (site 1)
% S1_SWS_total = timetable( ...
%     VWC1_mean.Time, ...
%     sum(VWC1_mean{:,layers},2,'omitnan')*fac, ... % total SWS (mean)
%     sum(VWC1_lo95{:,layers},2,'omitnan')*fac, ... % total SWS (min bound)
%     sum(VWC1_hi95{:,layers},2,'omitnan')*fac, ... % total SWS (max bound)
%     'VariableNames',{'SWS_total_mean_cm','SWS_total_min_cm','SWS_total_max_cm'});

S1_SWS_avg20 = timetable( ...
    VWC1_mean.Time, ...
    mean(VWC1_mean{:,layers},2,'omitnan')*fac, ... % avg per 20 cm (mean)
    mean(VWC1_lo95{:,layers},2,'omitnan')*fac, ... % avg per 20 cm (min bound)
    mean(VWC1_hi95{:,layers},2,'omitnan')*fac, ... % avg per 20 cm (max bound)
    'VariableNames',{'SWS_avg20_mean_cm','SWS_avg20_min_cm','SWS_avg20_max_cm'});

% ---------- Woodland (site 2)
% S2_SWS_total = timetable( ...
%     VWC2_mean.Time, ...
%     sum(VWC2_mean{:,layers},2,'omitnan')*fac, ...
%     sum(VWC2_lo95{:,layers},2,'omitnan')*fac, ...
%     sum(VWC2_hi95{:,layers},2,'omitnan')*fac, ...
%     'VariableNames',{'SWS_total_mean_cm','SWS_total_min_cm','SWS_total_max_cm'});

S2_SWS_avg20 = timetable( ...
    VWC2_mean.Time, ...
    mean(VWC2_mean{:,layers},2,'omitnan')*fac, ...
    mean(VWC2_lo95{:,layers},2,'omitnan')*fac, ...
    mean(VWC2_hi95{:,layers},2,'omitnan')*fac, ...
    'VariableNames',{'SWS_avg20_mean_cm','SWS_avg20_min_cm','SWS_avg20_max_cm'});

% above is SWS within 20 cm
%% average SWS every 20 cm -> top/mid/bot
% Requires: S1_SWS_avg20 and S2_SWS_avg20 from the previous step
avg1_mean = S1_SWS_avg20.SWS_avg20_mean_cm;
avg1_min  = S1_SWS_avg20.SWS_avg20_min_cm;
avg1_max  = S1_SWS_avg20.SWS_avg20_max_cm;

avg2_mean = S2_SWS_avg20.SWS_avg20_mean_cm;
avg2_min  = S2_SWS_avg20.SWS_avg20_min_cm;
avg2_max  = S2_SWS_avg20.SWS_avg20_max_cm;

% Depth-interval multipliers → total factors
S1_factor = 15*1.00 + 8*1.80 + 3*3.12;   % = 38.76   (0–520 cm)
S2_factor = 16*1.00 + 3*2.74 + 3*4.71;   % = 38.35   (0–440 cm)

% Totals
S1_SWS_total = timetable( ...
    S1_SWS_avg20.Time, ...
    avg1_mean*S1_factor, avg1_min*S1_factor, avg1_max*S1_factor, ...
    'VariableNames', {'SWS_0_520_total_mean_cm','SWS_0_520_total_min_cm','SWS_0_520_total_max_cm'});

S2_SWS_total = timetable( ...
    S2_SWS_avg20.Time, ...
    avg2_mean*S2_factor, avg2_min*S2_factor, avg2_max*S2_factor, ...
    'VariableNames', {'SWS_0_440_total_mean_cm','SWS_0_440_total_min_cm','SWS_0_440_total_max_cm'});

%% daily change in SWS

% Savanna: 0–520 cm
S1 = sortrows(S1_SWS_total);
S1_dSWS = timetable( ...
    S1.Time(2:end), ...
    diff(S1.SWS_0_520_total_mean_cm), ...
    diff(S1.SWS_0_520_total_min_cm), ...
    diff(S1.SWS_0_520_total_max_cm), ...
    'VariableNames', {'dSWS_mean','dSWS_min','dSWS_max'});

% Woodland: 0–440 cm
S2 = sortrows(S2_SWS_total);
S2_dSWS = timetable( ...
    S2.Time(2:end), ...
    diff(S2.SWS_0_440_total_mean_cm), ...
    diff(S2.SWS_0_440_total_min_cm), ...
    diff(S2.SWS_0_440_total_max_cm), ...
    'VariableNames', {'dSWS_mean','dSWS_min','dSWS_max'});

% (Optional) quick peek
head(S1_dSWS), head(S2_dSWS)

%% time series of std of ΔSWS
diff_S1mean_min = S1_dSWS.dSWS_mean-S1_dSWS.dSWS_min;
diff_S1mean_max = S1_dSWS.dSWS_max-S1_dSWS.dSWS_mean;
std_S1 = diff_S1mean_max/1.96;

diff_S2mean_max = S2_dSWS.dSWS_max-S2_dSWS.dSWS_mean;
std_S2 = diff_S2mean_max/1.96;

%% absolute uncertainty of ΔGW
% absolute uncertainty of 10-min rainfall rate from tipping bucket rain gauge
P_WB = daily_rain_filter(2:85,:);

% assume the uncertainty of 10-min P measured from TB RG = 0.1 (10%)
std_P_WB = 0.1*table2array(P_WB(:,2));

% absolute uncertainty of ET from openET -- RMSE in cm!!
ETgrass_WB = ET_P_grass2024(2:85,4);
ETwood_WB = ET_P_wood2024(2:85,4);

rmse_grassland = 1.08*0.1;
rmse_shrubland = 0.84*0.1;

% water balance for ΔGW = P-ET-dSWS -- S1:savanna; S2:woodland
dGW_S1 = table2array(P_WB(:,2)) - table2array(ETgrass_WB) - table2array(S1_dSWS(:,1));
dGW_S2 = table2array(P_WB(:,2)) - table2array(ETwood_WB) - table2array(S2_dSWS(:,1));

dGW_S1_raw = dGW_S1;
dGW_S2_raw = dGW_S2;


% optional smoothing / percolation delay
alpha = 1;   % 0 < alpha <= 1 ; lower = more lag, higher = closer to original
for k = 2:numel(dGW_S1)
    dGW_S1(k) = (1-alpha)*dGW_S1(k-1) + alpha*dGW_S1(k);
    dGW_S2(k) = (1-alpha)*dGW_S2(k-1) + alpha*dGW_S2(k);
end

% uncertainty of dGW
std_dGW_S1 = (std_P_WB.^2 + rmse_grassland.^2 + std_S1.^2).^0.5;
std_dGW_S2 = (std_P_WB.^2 + rmse_shrubland.^2 + std_S2.^2).^0.5;

% 95%CI of dGW
CI95_dGW_S1 = 1.96*std_dGW_S1;
CI95_dGW_S2 = 1.96*std_dGW_S2;

%% truc dGW
GW1 = GW_series{1,2};              % Savanna depth (cm)
GW2 = GW_series{1,4};              % Woodland depth (cm)

% Daily change with sign convention: positive = rise (depth decreases)
dGW_true_S1 = -(GW1(2:85) - GW1(1:84));   % length 84
dGW_true_S2 = -(GW2(2:85) - GW2(1:84));   % length 84

%% keep the direct ouput of GWL change i.e., q
% time
t  = P_WB{:,1}; 
tn = datenum(t);
tr = timerange(t(1), t(end), 'closed');

% q (method 4) aligned
q_s1 = q_combined_site1(tr,'q_method4').q_method4;
q_s2 = q_combined_site2(tr,'q_method4').q_method4;

% components (cm d^-1)
P_cm   = P_WB.Rain_cm_Tot;
ETg_cm = ETgrass_WB.ET_cm;
ETw_cm = ETwood_WB.ET_cm;
dSWS1  = S1_dSWS.dSWS_mean;
dSWS2  = S2_dSWS.dSWS_mean;

% convert ET to mm for left axis in components panels
ETg_mm = ETg_cm*10;
ETw_mm = ETw_cm*10;

% --- uncertainties (must be 84x1, aligned to t) ---
% P in cm d^-1
sigmaP_cm = std_P_WB;                      % already cm/day, 84x1

% ET RMSE: defaults if not already defined (OpenET daily)
if ~exist('rmse_grassland','var'),  rmse_grassland  = 1.08/10; end  % unit:cm/d
if ~exist('rmse_shrubland','var'),  rmse_shrubland  = 0.84/10; end  % unit:cm/d
sigmaETg_mm = (rmse_grassland*10) * ones(size(ETg_cm));   % mm/d
sigmaETw_mm = (rmse_shrubland*10) * ones(size(ETw_cm));   % mm/d

% ΔSWS in cm d^-1
sigmaS1_cm = std_S1;                        % cm/d, 84x1
sigmaS2_cm = std_S2;                        % cm/d, 84x1


% colors
colWB   = [0 0.447 0.741];
colMeas = [0.85 0.1 0.1];
colQ    = [0.15 0.15 0.15];
colPbar = [0.2 0.6 0.9];
colET   = [0.20 0.60 0.20];
colSWS  = [0.85 0.33 0.10];

figure('Position',[80 80 1150 700]);
tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

% ---------- A) Savanna: ΔGW & q ----------
ax1 = nexttile(3); hold on
fill([tn; flipud(tn)], [dGW_S1-CI95_dGW_S1; flipud(dGW_S1+CI95_dGW_S1)], ...
     colWB, 'FaceAlpha',0.15, 'EdgeColor','none');
plot(tn, dGW_S1, 'LineWidth',1.6, 'Color',colWB);
plot(tn, dGW_true_S1, 'Color',colMeas, 'LineWidth',1.5);
plot(tn, q_s1, '--', 'Color',colQ, 'LineWidth',1.2);
yline(0,'k:');
ylabel('\DeltaGWL & drainage (cm d^{-1})');xlabel('Date');
title('Savanna');
legend('WB based \DeltaGWL 95% CI','WB based \DeltaGWL','Measured \DeltaGWL','Drainage q', ...
       'Location','best');

% ---------- B) Woodland: ΔGW & q ----------
ax2 = nexttile(4); hold on
fill([tn; flipud(tn)], [dGW_S2-CI95_dGW_S2; flipud(dGW_S2+CI95_dGW_S2)], ...
     colWB, 'FaceAlpha',0.15, 'EdgeColor','none');
plot(tn, dGW_S2, 'LineWidth',1.6, 'Color',colWB);
plot(tn, dGW_true_S2, 'Color',colMeas, 'LineWidth',1.5);
plot(tn, q_s2, '--', 'Color',colQ, 'LineWidth',1.2);
yline(0,'k:');
ylabel('\DeltaGWL & drainage (cm d^{-1})');xlabel('Date');
title('Woodland');

% --- make the top two panels share identical y-limits ---
yl_all = [dGW_S1-CI95_dGW_S1; dGW_S1+CI95_dGW_S1; dGW_true_S1; ...
          dGW_S2-CI95_dGW_S2; dGW_S2+CI95_dGW_S2; dGW_true_S2];
ylmin = min(yl_all); ylmax = max(yl_all);
pad   = 0.05*(ylmax-ylmin + eps);
set([ax1 ax2],'YLim',[ylmin-pad ylmax+pad]);

% ---------- C) Savanna components (dual y-axes) ----------
ax3 = nexttile(1); hold on   % was: ax3 = nexttile;
yyaxis left
b1 = bar(tn, P_cm, 'FaceColor',colPbar, 'FaceAlpha',0.5, ...
         'EdgeColor','none','BarWidth',1); uistack(b1,'bottom');
% 'Color',[0 0.45 0.74]
eP  = errorbar(tn, P_cm, sigmaP_cm, 'LineStyle','none', ...
               'Color',[0 0.45 0.74], 'CapSize',0, 'LineWidth',1, ...
               'HandleVisibility','off');              % << hide in legend
ETg_hi = ETg_mm + sigmaETg_mm;  ETg_lo = ETg_mm - sigmaETg_mm;
fill([tn; flipud(tn)], [ETg_lo; flipud(ETg_hi)], colET, ...
     'FaceAlpha',0.15, 'EdgeColor','none', 'HandleVisibility','off');
hETg = plot(tn, ETg_mm, 'Color',colET, 'LineWidth',1.5,'LineStyle', '--');
ylabel('P (cm d^{-1}), ET (mm d^{-1})');

yyaxis right
S1_hi = dSWS1 + sigmaS1_cm;  S1_lo = dSWS1 - sigmaS1_cm;
fill([tn; flipud(tn)], [S1_lo; flipud(S1_hi)], colSWS, ...
     'FaceAlpha',0.15, 'EdgeColor','none', 'HandleVisibility','off');
hSWS1 = plot(tn, dSWS1, 'Color',colSWS, 'LineWidth',1.5,'LineStyle', '-');
ylabel('\DeltaSWS (cm d^{-1})');
title('Savanna');
legend([b1 hETg hSWS1], {'P','ET','\DeltaSWS'}, 'Location','best');  % << only these

% ---------- D) Woodland components (dual y-axes) ----------
ax4 = nexttile(2); hold on
yyaxis left
b2 = bar(tn, P_cm, 'FaceColor',colPbar, 'FaceAlpha',0.5, ...
         'EdgeColor','none','BarWidth',1); uistack(b2,'bottom');
ePw = errorbar(tn, P_cm, sigmaP_cm, 'LineStyle','none', ...
               'Color',[0 0.45 0.74], 'CapSize',0, 'LineWidth',1, ...
               'HandleVisibility','off');              % << hide in legend
ETw_hi = ETw_mm + sigmaETw_mm;  ETw_lo = ETw_mm - sigmaETw_mm;
fill([tn; flipud(tn)], [ETw_lo; flipud(ETw_hi)], colET, ...
     'FaceAlpha',0.15, 'EdgeColor','none', 'HandleVisibility','off');
hETw = plot(tn, ETw_mm, 'Color',colET, 'LineWidth',1.5,'LineStyle', '--');
ylabel('P (cm d^{-1}), ET (mm d^{-1})'); 

yyaxis right
S2_hi = dSWS2 + sigmaS2_cm;  S2_lo = dSWS2 - sigmaS2_cm;
fill([tn; flipud(tn)], [S2_lo; flipud(S2_hi)], colSWS, ...
     'FaceAlpha',0.15, 'EdgeColor','none', 'HandleVisibility','off');
hSWS2 = plot(tn, dSWS2, 'Color',colSWS, 'LineWidth',1.5,'LineStyle', '-');
ylabel('\DeltaSWS (cm d^{-1})');
title('Woodland');
% legend([b2 hETw hSWS2], {'P','ET','\DeltaSWS'}, 'Location','best');  % << only these


% --- Make the right y-axes (ΔSWS) identical on the bottom two panels ---
sws_all = [dSWS1(:)-sigmaS1_cm(:); dSWS1(:)+sigmaS1_cm(:); ...
           dSWS2(:)-sigmaS2_cm(:); dSWS2(:)+sigmaS2_cm(:)];
rmin = min(sws_all); rmax = max(sws_all);
pad  = 0.05*(rmax - rmin + eps);

yyaxis(ax3,'right'); ylim(ax3, [rmin - pad, rmax + pad]);
yyaxis(ax4,'right'); ylim(ax4, [rmin - pad, rmax + pad]);
 

% ----- LEFT y-axis (P in cm d^-1, ET in mm d^-1): start at 0, shared top -----
Lmax = max([ P_cm(:) + sigmaP_cm(:); ...
             ETg_mm(:) + sigmaETg_mm(:); ...
             ETw_mm(:) + sigmaETw_mm(:) ]);
Lpad = 0.05*(Lmax + eps);                 % small headroom

yyaxis(ax3,'left');  ylim(ax3, [0, Lmax + Lpad]);
yyaxis(ax4,'left');  ylim(ax4, [0, Lmax + Lpad]);


% ---- consistent x-limits / date ticks on all panels ----
axs = [ax1 ax2 ax3 ax4];
tick_dt  = t(1):caldays(14):t(end);        % weekly ticks (change if desired)
tick_num = datenum(tick_dt);
for ax = axs
    set(ax, 'XLim', [tn(1) tn(end)], 'XTick', tick_num);
    datetick(ax, 'x', 'mmm dd', 'keepticks','keeplimits'); % keep date strings
    ax.XAxis.Exponent = 0;
    grid(ax,'on'); box(ax,'on');
end

linkaxes(axs,'x');   % synchronized zoom/pan

% After you've finished drawing ax3 and ax4 (the components panels)
for ax = [ax3 ax4]
    % left y-axis (P & ET)
    yyaxis(ax,'left');
    ax.YAxis(1).Color = 'k';            % tick/axis color
    ax.YAxis(1).Label.Color = 'k';      % label color

    % right y-axis (ΔSWS)
    yyaxis(ax,'right');
    ax.YAxis(2).Color = 'k';
    ax.YAxis(2).Label.Color = 'k';
end

%% on how many days measured and estimated drainage mostly fell within the 95% CI of water balance-based ΔGWL
% Helper to count "within CI" days
withinCI = @(x,mu,ci) x >= (mu - ci) & x <= (mu + ci);

% ---------- Savanna (site 1) ----------
lb1 = dGW_S1 - CI95_dGW_S1;
ub1 = dGW_S1 + CI95_dGW_S1;

% Measured ΔGWL vs WB CI
valid_meas1 = ~isnan(dGW_true_S1) & ~isnan(lb1) & ~isnan(ub1);
in_meas1    = withinCI(dGW_true_S1, dGW_S1, CI95_dGW_S1) & valid_meas1;
n_meas1     = sum(in_meas1);
N_meas1     = sum(valid_meas1);
pct_meas1   = 100*n_meas1/N_meas1;

% 1-m drainage q vs WB CI
valid_q1  = ~isnan(q_s1(:,1)) & ~isnan(lb1) & ~isnan(ub1);
in_q1     = withinCI(q_s1(:,1), dGW_S1, CI95_dGW_S1) & valid_q1;
n_q1      = sum(in_q1);
N_q1      = sum(valid_q1);
pct_q1    = 100*n_q1/N_q1;

fprintf('Savanna: measured ΔGWL within CI on %d/%d days (%.1f%%)\n', n_meas1, N_meas1, pct_meas1);
fprintf('Savanna: drainage q within CI on %d/%d days (%.1f%%)\n',   n_q1,    N_q1,    pct_q1);

% ---------- Woodland (site 2) ----------
lb2 = dGW_S2 - CI95_dGW_S2;
ub2 = dGW_S2 + CI95_dGW_S2;

% Measured ΔGWL vs WB CI
valid_meas2 = ~isnan(dGW_true_S2) & ~isnan(lb2) & ~isnan(ub2);
in_meas2    = withinCI(dGW_true_S2, dGW_S2, CI95_dGW_S2) & valid_meas2;
n_meas2     = sum(in_meas2);
N_meas2     = sum(valid_meas2);
pct_meas2   = 100*n_meas2/N_meas2;

% 1-m drainage q vs WB CI
valid_q2  = ~isnan(q_s2(:,1)) & ~isnan(lb2) & ~isnan(ub2);
in_q2     = withinCI(q_s2(:,1), dGW_S2, CI95_dGW_S2) & valid_q2;
n_q2      = sum(in_q2);
N_q2      = sum(valid_q2);
pct_q2    = 100*n_q2/N_q2;

fprintf('Woodland: measured ΔGWL within CI on %d/%d days (%.1f%%)\n', n_meas2, N_meas2, pct_meas2);
fprintf('Woodland: drainage q within CI on %d/%d days (%.1f%%)\n',   n_q2,    N_q2,    pct_q2);

%% table summary: table 5
% P
CI95_P = 1.96*sigmaP_cm; % cm
CI95_P_up = P_cm + CI95_P;
CI95_P_lo = P_cm - CI95_P;

P_mean = nanmean(P_cm)
P_mean_up = nanmean(CI95_P_up)
P_mean_lo = nanmean(CI95_P_lo)

% ET
CI95_ETg = 1.96*rmse_grassland; % cm
CI95_ETw = 1.96*rmse_shrubland; % cm
 
CI95_ETg_up = ETg_cm + CI95_ETg;
CI95_ETg_lo = ETg_cm - CI95_ETg;

CI95_ETw_up = ETw_cm + CI95_ETw;
CI95_ETw_lo = ETw_cm - CI95_ETw;

ETg_mean = nanmean(ETg_cm)
ETg_mean_up = nanmean(CI95_ETg_up)
ETg_mean_lo = nanmean(CI95_ETg_lo)

ETw_mean = nanmean(ETw_cm)
ETw_mean_up = nanmean(CI95_ETw_up)
ETw_mean_lo = nanmean(CI95_ETw_lo)

% SWS

SWSg_mean = nanmean(S1_dSWS.dSWS_mean)
SWSg_mean_up = nanmean(S1_dSWS.dSWS_max)
SWSg_mean_lo = nanmean(S1_dSWS.dSWS_min)

SWSw_mean = nanmean(S2_dSWS.dSWS_mean)
SWSw_mean_up = nanmean(S2_dSWS.dSWS_max)
SWSw_mean_lo = nanmean(S2_dSWS.dSWS_min)

% GWL

dGW_S1_mean = nanmean(dGW_S1)
dGW_S1_mean_up = nanmean(dGW_S1+CI95_dGW_S1)
dGW_S1_mean_lo = nanmean(dGW_S1-CI95_dGW_S1)

dGW_S2_mean = nanmean(dGW_S2)
dGW_S2_mean_up = nanmean(dGW_S2+CI95_dGW_S2)
dGW_S2_mean_lo = nanmean(dGW_S2-CI95_dGW_S2)

