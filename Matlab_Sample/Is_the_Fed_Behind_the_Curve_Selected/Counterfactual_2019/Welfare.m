cd("..")
load("Welfare_average.mat");

temp = reshape(sum(Ey_gap_raw.^2, 2), [], 3);
Loss_raw = sum(temp(:, 1:2), 2);

temp = reshape(sum(Ey_gap.^2, 2), [], 3);
Loss = sum(temp(:, 1:2), 2);

% Calculate Welfare_loss%
Welfare_loss =  100*(Loss_raw - Loss)./Loss;

Welfare_accum = zeros(length(Loss), 1);

for i = 1:length(Loss(:, 1))
    Welfare_accum(i:length(Loss(:, 1)), 1) = Welfare_accum(i:length(Loss(:, 1)), 1)+(Loss_raw(i, 1) - Loss(i, 1));
end

Welfare_accum_perc = zeros(length(Loss), 1);

for i = 1:length(Loss(:, 1))
    Welfare_accum_perc(i:length(Loss(:, 1)), 1) = Welfare_accum_perc(i:length(Loss(:, 1)), 1)+(Loss_raw(i, 1) - Loss(i, 1))./Loss(i, 1);
end


%% FIGURE
figure;
hold on;

% Set gradient color
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
xlim([2019, 2024.25]);
xlabel('Year');
ylabel('Welfare improved\%');
hold off;