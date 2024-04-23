%% Baseline Welfare Improved
cd("Counterfactual_2019/")
load("Welfare_foreadj2.mat")

temp = reshape(sum(Ey_gap_raw.^2, 2), [], 3);
Loss_raw = sum(temp(:, 1:2), 2);

temp = reshape(sum(Ey_gap.^2, 2), [], 3);
Loss = sum(temp(:, 1:2), 2);

% Calculate Welfare_loss%
Welfare_loss =  100*(Loss_raw - Loss)./Loss;


%% Welfare Improved H20

load("Welfare_Hf20.mat")

temp = reshape(sum(Ey_gap_raw.^2, 2), [], 3);
Loss_raw_Hf20 = sum(temp(:, 1:2), 2);

temp = reshape(sum(Ey_gap.^2, 2), [], 3);
Loss_Hf20 = sum(temp(:, 1:2), 2);

% Calculate Welfare_loss%
Welfare_loss_Hf20 =  100*(Loss_raw_Hf20 - Loss_Hf20)./Loss_Hf20;

%% Welfare Improved H12

load("Welfare_Hf12.mat")

temp = reshape(sum(Ey_gap_raw.^2, 2), [], 3);
Loss_raw_Hf12 = sum(temp(:, 1:2), 2);

temp = reshape(sum(Ey_gap.^2, 2), [], 3);
Loss_Hf12 = sum(temp(:, 1:2), 2);

% Calculate Welfare_loss%
Welfare_loss_Hf12 =  100*(Loss_raw_Hf12 - Loss_Hf12)./Loss_Hf12;

%% Welfare Improved H8

load("Welfare_Hf8.mat")

temp = reshape(sum(Ey_gap_raw.^2, 2), [], 3);
Loss_raw_Hf8 = sum(temp(:, 1:2), 2);

temp = reshape(sum(Ey_gap.^2, 2), [], 3);
Loss_Hf8 = sum(temp(:, 1:2), 2);

% Calculate Welfare_loss%
Welfare_loss_Hf8 =  100*(Loss_raw_Hf8 - Loss_Hf8)./Loss_Hf8;

%% Welfare Improved H0

load("Welfare_Hf0.mat")

temp = reshape(sum(Ey_gap_raw.^2, 2), [], 3);
Loss_raw_Hf0 = sum(temp(:, 1:2), 2);

temp = reshape(sum(Ey_gap.^2, 2), [], 3);
Loss_Hf0 = sum(temp(:, 1:2), 2);

% Calculate Welfare_loss%
Welfare_loss_Hf0 =  100*(Loss_raw_Hf0 - Loss_Hf0)./Loss_Hf0;

%% Figure Baseline

figure;
subplot(1, 5, 1);
hold on;

% Gradient color
startColor_p = [0.7, 0.7, 1]; % Positive value, light blue
endColor_p = [0, 0, 0.5]; % Positive value, dark blue

startColor_n = [1, 0.7, 0.7]; % Negative value, light red
endColor_n = [0.5, 0, 0]; % Negative value, dark red

for i = 1:numel(Welfare_loss)
    % Calculate the color of the bar
    if Welfare_loss(i) >= 0
        fraction = abs(Welfare_loss(i)) / max(abs(Welfare_loss));
        barColor = startColor_p + (endColor_p - startColor_p) * fraction;

        bar(yearn(i), Welfare_loss(i), 'FaceColor', barColor, 'BarWidth', 0.2, 'EdgeColor', 'none');
    else
        fraction = abs(Welfare_loss(i)) / max(abs(Welfare_loss));
        barColor = startColor_n + (endColor_n - startColor_n) * fraction;

        bar(yearn(i), Welfare_loss(i), 'FaceColor', barColor, 'BarWidth', 0.2, 'EdgeColor', 'none');
    end
end

set(0,'defaulttextInterpreter','latex') 
xlim([2019, 2024]);
ylim([-50, 110]);
xlabel('Year');
ylabel('Welfare improved\%');
hold off

%% Figure Hf20

subplot(1, 5, 2);
hold on;

% Gradient color
startColor_p = [0.7, 0.7, 1]; % Positive value, light blue
endColor_p = [0, 0, 0.5]; % Positive value, dark blue

startColor_n = [1, 0.7, 0.7]; % Negative value, light red
endColor_n = [0.5, 0, 0]; % Negative value, dark red

for i = 1:numel(Welfare_loss_Hf20)
    % Calculate the color of the bar
    if Welfare_loss_Hf20(i) >= 0
        fraction = abs(Welfare_loss_Hf20(i)) / max(abs(Welfare_loss_Hf20));
        barColor = startColor_p + (endColor_p - startColor_p) * fraction;

        bar(yearn(i), Welfare_loss_Hf20(i), 'FaceColor', barColor, 'BarWidth', 0.2, 'EdgeColor', 'none');
    else
        fraction = abs(Welfare_loss_Hf20(i)) / max(abs(Welfare_loss_Hf20));
        barColor = startColor_n + (endColor_n - startColor_n) * fraction;

        bar(yearn(i), Welfare_loss_Hf20(i), 'FaceColor', barColor, 'BarWidth', 0.2, 'EdgeColor', 'none');
    end
end

set(0,'defaulttextInterpreter','latex') 
xlim([2019, 2024]);
ylim([-50, 110]);
xlabel('Year');
ylabel('Welfare improved\% Hf20');
hold off


%% Figure Hf12

subplot(1, 5, 3);
hold on;

% Gradient color
startColor_p = [0.7, 0.7, 1]; % Positive value, light blue
endColor_p = [0, 0, 0.5]; % Positive value, dark blue

startColor_n = [1, 0.7, 0.7]; % Negative value, light red
endColor_n = [0.5, 0, 0]; % Negative value, dark red

for i = 1:numel(Welfare_loss_Hf12)
    % Calculate the color of the bar
    if Welfare_loss_Hf12(i) >= 0
        fraction = abs(Welfare_loss_Hf12(i)) / max(abs(Welfare_loss_Hf12));
        barColor = startColor_p + (endColor_p - startColor_p) * fraction;

        bar(yearn(i), Welfare_loss_Hf12(i), 'FaceColor', barColor, 'BarWidth', 0.2, 'EdgeColor', 'none');
    else
        fraction = abs(Welfare_loss_Hf12(i)) / max(abs(Welfare_loss_Hf12));
        barColor = startColor_n + (endColor_n - startColor_n) * fraction;

        bar(yearn(i), Welfare_loss_Hf12(i), 'FaceColor', barColor, 'BarWidth', 0.2, 'EdgeColor', 'none');
    end
end

set(0,'defaulttextInterpreter','latex') 
xlim([2019, 2024]);
ylim([-50, 110]);
xlabel('Year');
ylabel('Welfare improved\% Hf12');
hold off


%% Figure Hf8

subplot(1, 5, 4);
hold on;

% Gradient color
startColor_p = [0.7, 0.7, 1]; % Positive value, light blue
endColor_p = [0, 0, 0.5]; % Positive value, dark blue

startColor_n = [1, 0.7, 0.7]; % Negative value, light red
endColor_n = [0.5, 0, 0]; % Negative value, dark red

for i = 1:numel(Welfare_loss_Hf8)
    % Calculate the color of the bar
    if Welfare_loss_Hf8(i) >= 0
        fraction = abs(Welfare_loss_Hf8(i)) / max(abs(Welfare_loss_Hf8));
        barColor = startColor_p + (endColor_p - startColor_p) * fraction;

        bar(yearn(i), Welfare_loss_Hf8(i), 'FaceColor', barColor, 'BarWidth', 0.2, 'EdgeColor', 'none');
    else
        fraction = abs(Welfare_loss_Hf8(i)) / max(abs(Welfare_loss_Hf8));
        barColor = startColor_n + (endColor_n - startColor_n) * fraction;

        bar(yearn(i), Welfare_loss_Hf8(i), 'FaceColor', barColor, 'BarWidth', 0.2, 'EdgeColor', 'none');
    end
end

set(0,'defaulttextInterpreter','latex') 
xlim([2019, 2024]);
ylim([-50, 110]);
xlabel('Year');
ylabel('Welfare improved\% Hf8');
hold off


%% Figure Hf0

subplot(1, 5, 5);
hold on;

% Gradient color
startColor_p = [0.7, 0.7, 1]; % Positive value, light blue
endColor_p = [0, 0, 0.5]; % Positive value, dark blue

startColor_n = [1, 0.7, 0.7]; % Negative value, light red
endColor_n = [0.5, 0, 0]; % Negative value, dark red

for i = 1:numel(Welfare_loss_Hf0)
    % Calculate the color of the bar
    if Welfare_loss_Hf0(i) >= 0
        fraction = abs(Welfare_loss_Hf0(i)) / max(abs(Welfare_loss_Hf0));
        barColor = startColor_p + (endColor_p - startColor_p) * fraction;

        bar(yearn(i), Welfare_loss_Hf0(i), 'FaceColor', barColor, 'BarWidth', 0.2, 'EdgeColor', 'none');
    else
        fraction = abs(Welfare_loss_Hf0(i)) / max(abs(Welfare_loss_Hf0));
        barColor = startColor_n + (endColor_n - startColor_n) * fraction;

        bar(yearn(i), Welfare_loss_Hf0(i), 'FaceColor', barColor, 'BarWidth', 0.2, 'EdgeColor', 'none');
    end
end

set(0,'defaulttextInterpreter','latex') 
xlim([2019, 2024]);
ylim([-50, 180]);
xlabel('Year');
ylabel('Welfare improved\% No adjustment');
hold off