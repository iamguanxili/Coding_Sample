load("bTable_output_SEP_2019.mat")

% Median and confidence intervals of counterfactual paths
UR_output = squeeze(bTable_output_SEP(:, 4, :));
PI_output = squeeze(bTable_output_SEP(:, 5, :));
FFR_output = squeeze(bTable_output_SEP(:, 6, :));

UR_output_lb1 = prctile(UR_output, 10, 2);
UR_output_ub1 = prctile(UR_output, 100-10, 2);
UR_output_lb2 = prctile(UR_output, 20, 2);
UR_output_ub2 = prctile(UR_output, 100-20, 2);
UR_output_med = prctile(UR_output, 50, 2);

PI_output_lb1 = prctile(PI_output, 10, 2);
PI_output_ub1 = prctile(PI_output, 100-10, 2);
PI_output_lb2 = prctile(PI_output, 20, 2);
PI_output_ub2 = prctile(PI_output, 100-20, 2);
PI_output_med = prctile(PI_output, 50, 2);

FFR_output_lb1 = prctile(FFR_output, 10, 2);
FFR_output_ub1 = prctile(FFR_output, 100-10, 2);
FFR_output_lb2 = prctile(FFR_output, 20, 2);
FFR_output_ub2 = prctile(FFR_output, 100-20, 2);
FFR_output_med = prctile(FFR_output, 50, 2);

% Long-run targets
UR_raw_table = readtable('fomc_projections_MYrep.xlsx', 'Sheet', 2);
PI_raw_table = readtable('fomc_projections_MYrep.xlsx', 'Sheet', 1);

UR_raw_table.median_lr(isnan(UR_raw_table.median_lr)) = UR_raw_table.RealTime_LR(isnan(UR_raw_table.median_lr));
PI_raw_table.median_lr(isnan(PI_raw_table.median_lr)) = PI_raw_table.RealTime_LR(isnan(PI_raw_table.median_lr));

UR_raw_table = UR_raw_table;
PI_raw_table = PI_raw_table;


% Raw data
rawbTable = readtable('Counterfactual_Objectives.xlsx');

yearn_raw = rawbTable.yearn;
UR_raw = rawbTable.UR;
PI_raw = rawbTable.PI;
FFR_raw = rawbTable.FFR;

% Create tables to store counterfactuals
% UR
URtable_output = table(yearn_raw);

URtable_output.UR_output_lb1 = UR_output_lb1;
URtable_output.UR_output_ub1 = UR_output_ub1;
URtable_output.UR_output_lb2 = UR_output_lb2;
URtable_output.UR_output_ub2 = UR_output_ub2;
URtable_output.UR_output_med = UR_output_med;

% PI
PItable_output = table(yearn_raw); % 使用 rawTable 中的 yearn 列作为第一列

PItable_output.PI_output_lb1 = PI_output_lb1;
PItable_output.PI_output_ub1 = PI_output_ub1;
PItable_output.PI_output_lb2 = PI_output_lb2;
PItable_output.PI_output_ub2 = PI_output_ub2;
PItable_output.PI_output_med = PI_output_med;

% FFR
FFRtable_output = table(yearn_raw); % 使用 rawTable 中的 yearn 列作为第一列

FFRtable_output.FFR_output_lb1 = FFR_output_lb1;
FFRtable_output.FFR_output_ub1 = FFR_output_ub1;
FFRtable_output.FFR_output_lb2 = FFR_output_lb2;
FFRtable_output.FFR_output_ub2 = FFR_output_ub2;
FFRtable_output.FFR_output_med = FFR_output_med;

%% Figure

% PI
figure;
set(0,'defaulttextInterpreter','latex')
subplot(2,2,1)
hold on;

zlb_init = min(find(yearn >= 2008.75));
zlb_end = min(find(yearn >= 2015.6));
zlb2_init = min(find(yearn >= 2020.2));
zlb2_end = min(find(yearn >= 2021.8));

rectangle('Position', ...
[yearn(zlb_init) -30 yearn(zlb_end)-yearn(zlb_init) 100], ...
'FaceColor', [0.957, 0.643, 0.176, 0.1], 'EdgeColor', 'none')
rectangle('Position', ...
[yearn(zlb2_init) -30 yearn(zlb2_end)-yearn(zlb2_init) 100], ...
'FaceColor', [0.957, 0.643, 0.176, 0.1], 'EdgeColor', 'none')

plot(yearn_raw, PI_raw, 'r-.', 'LineWidth', 1.5, 'DisplayName', 'raw\_PI');
plot(yearn_raw, PItable_output.PI_output_med, 'k', 'LineWidth', 1.5, 'DisplayName', 'adjusted\_PI');
plot(PI_raw_table.year_1, PI_raw_table.median_lr, 'Color', [0.7, 0.7, 0.7], 'LineStyle', '-.', 'DisplayName', 'longrun\_target');
u = fill([yearn_raw; flipud(yearn_raw)], [PItable_output.PI_output_lb1; flipud(PItable_output.PI_output_ub1)], [0.8, 0.8, 0.8], 'LineStyle', 'none');
l = fill([yearn_raw; flipud(yearn_raw)], [PItable_output.PI_output_lb2; flipud(PItable_output.PI_output_ub2)], [0.4, 0.4, 0.4], 'LineStyle', 'none');
set(u, 'HandleVisibility', 'off');
set(l, 'HandleVisibility', 'off');
u.FaceAlpha = 0.3;
l.FaceAlpha = 0.1;
xlabel('Year');
ylabel('PI');
legend('show');
ylim([-2, 8]);
xlim([2020.75, 2023.75]);
hold off;


% UR
subplot(2,2,2)
hold on;

zlb_init = min(find(yearn >= 2008.75));
zlb_end = min(find(yearn >= 2015.6));
zlb2_init = min(find(yearn >= 2020.2));
zlb2_end = min(find(yearn >= 2021.8));

rectangle('Position', ...
[yearn(zlb_init) -30 yearn(zlb_end)-yearn(zlb_init) 100], ...
'FaceColor', [0.957, 0.643, 0.176, 0.1], 'EdgeColor', 'none')
rectangle('Position', ...
[yearn(zlb2_init) -30 yearn(zlb2_end)-yearn(zlb2_init) 100], ...
'FaceColor', [0.957, 0.643, 0.176, 0.1], 'EdgeColor', 'none')
plot(yearn_raw, UR_raw, 'r-.', 'LineWidth', 1.5, 'DisplayName', 'raw\_UR');
plot(yearn_raw, URtable_output.UR_output_med, 'k', 'LineWidth', 1.5, 'DisplayName', 'adjusted\_UR');
plot(UR_raw_table.year_1, UR_raw_table.median_lr, 'Color', [0.7, 0.7, 0.7], 'LineStyle', '-.', 'DisplayName', 'longrun\_target');
u = fill([yearn_raw; flipud(yearn_raw)], [URtable_output.UR_output_lb1; flipud(URtable_output.UR_output_ub1)], [0.8, 0.8, 0.8], 'LineStyle', 'none');
l = fill([yearn_raw; flipud(yearn_raw)], [URtable_output.UR_output_lb2; flipud(URtable_output.UR_output_ub2)], [0.4, 0.4, 0.4], 'LineStyle', 'none');
set(u, 'HandleVisibility', 'off');
set(l, 'HandleVisibility', 'off');
u.FaceAlpha = 0.3;
l.FaceAlpha = 0.1;
xlabel('Year');
ylabel('UR');
legend('show');
xlim([2020.75, 2023.75]);
ylim([-3, 20]);
hold off;


% FFR
subplot(2,2,3);
hold on;

zlb_init = min(find(yearn >= 2008.75));
zlb_end = min(find(yearn >= 2015.6));
zlb2_init = min(find(yearn >= 2020.2));
zlb2_end = min(find(yearn >= 2021.8));

rectangle('Position', ...
[yearn(zlb_init) -30 yearn(zlb_end)-yearn(zlb_init) 100], ...
'FaceColor', [0.957, 0.643, 0.176, 0.1], 'EdgeColor', 'none')
rectangle('Position', ...
[yearn(zlb2_init) -30 yearn(zlb2_end)-yearn(zlb2_init) 100], ...
'FaceColor', [0.957, 0.643, 0.176, 0.1], 'EdgeColor', 'none')


plot(yearn_raw, FFR_raw, 'r-.', 'LineWidth', 1.5, 'DisplayName', 'raw\_FFR');
plot(yearn_raw, FFRtable_output.FFR_output_med, 'k', 'LineWidth', 1.5, 'DisplayName', 'adjusted\_FFR');
u = fill([yearn_raw; flipud(yearn_raw)], [FFRtable_output.FFR_output_lb1; flipud(FFRtable_output.FFR_output_ub1)], [0.8, 0.8, 0.8], 'LineStyle', 'none');
l = fill([yearn_raw; flipud(yearn_raw)], [FFRtable_output.FFR_output_lb2; flipud(FFRtable_output.FFR_output_ub2)], [0.4, 0.4, 0.4], 'LineStyle', 'none');
set(u, 'HandleVisibility', 'off');
set(l, 'HandleVisibility', 'off');
u.FaceAlpha = 0.3;
l.FaceAlpha = 0.1;
xlabel('Year');
ylabel('FFR');
legend('show');
xlim([2020.75, 2023.75]);
ylim([-3, 10]);
hold off;

%% Calculate deviations

PI_dev = PI_output_med - PI_raw;
UR_dev = UR_output_med - UR_raw;
FFR_dev = FFR_output_med - FFR_raw;

% Select data after 2019
FFR_dev_selected = FFR_dev(yearn_raw >= 2019);
PI_dev_selected = PI_dev(yearn_raw >= 2019);
UR_dev_selected = UR_dev(yearn_raw >= 2019);

year = yearn_raw(yearn_raw >= 2019);

% Create the table to store the results
T = table(year, FFR_dev_selected, PI_dev_selected, UR_dev_selected);

data = table2array(T);

data = round(data, 2);

fid = fopen('latex_table.tex', 'w');
fprintf(fid, '\\begin{table}[htbp]\n');
fprintf(fid, '\\centering\n');
fprintf(fid, '\\caption{Your caption here}\n');
fprintf(fid, '\\label{tab:my_label}\n');
fprintf(fid, '\\begin{tabular}{cccc}\n');
fprintf(fid, '\\toprule\n');
fprintf(fid, 'Year & FFR\_adj & PI\_adj & UR\_adj \\\\\n');
fprintf(fid, '\\midrule\n');
for i = 1:size(data, 1)
    fprintf(fid, '%.2f & %.2f & %.2f & %.2f \\\\\n', data(i, :));
end
fprintf(fid, '\\bottomrule\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\end{table}\n');
fclose(fid);