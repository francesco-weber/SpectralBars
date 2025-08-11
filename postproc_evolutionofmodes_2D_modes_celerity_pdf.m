clear; clc;
close all;

% Read intermediate results
t = load("res_interm_t.txt");
M = load("res_interm.txt");

% Set font size for all text elements
font_size = 12;  % You can adjust this value

% Initialize a variable to store celerity values at the last timestep
celerity_last_timestep = [];

figure
title('Temporal evolution of selected modes', 'Interpreter', 'latex', 'FontSize', font_size);

for i=7
    for j=1
        % Select the modes
        pp = (M(:,3)==i & M(:,4)==j);
        Xre = M(pp,1);
        Xim = M(pp,2);
        phi = atan2(Xim, Xre);
        amp = sqrt(Xre.^2 + Xim.^2);
        
        % Apply phase unwrapping to prevent discontinuities
        phi_unwrapped = unwrap(phi);
        
        % Compute the derivative of the unwrapped phase to get celerity
        celerity = [0; diff(phi_unwrapped)./diff(t(1:end,1))];

        % Apply amplitude threshold
        if amp<5e-3
            celerity=nan;
            phi=nan;
        end
                    

        % Store the celerity at the last time step
        if ~isempty(celerity)
            celerity_last_timestep = [celerity_last_timestep; celerity(end)/i];
        end

                                                                    
        % Plot amplitude
        subplot(3,1,1)
        plot(t(1:end,1), amp, 'DisplayName', ['Mode ', int2str(i), ',', int2str(j)]);
        ylabel('$a_{ij}$', 'Interpreter', 'latex', 'FontSize', font_size);
        %xlim([0 9000])
        hold on

        % Plot celerity, normalized by mode number i
        subplot(3,1,2)
        plot(t(1:end,1), celerity/i, 'DisplayName', ['Mode ', int2str(i), ',', int2str(j)]);
        ylabel('$\omega/\lambda$', 'Interpreter', 'latex', 'FontSize', font_size);
        xlabel('$\hat{t}$', 'Interpreter', 'latex', 'FontSize', font_size);
        %xlim([0 9000])
        hold on
    end
end

% Plot the distribution of celerity values at the last timestep
subplot(3,1,3)
histogram(celerity_last_timestep, 'Normalization', 'pdf','NumBins',10);
ylabel('pdf', 'Interpreter', 'latex', 'FontSize', font_size);
xlabel('$\omega/\lambda$', 'Interpreter', 'latex', 'FontSize', font_size);

% Set x and y ticks for all subplots
for k=1:3
    subplot(3,1,k)
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', font_size); % Use LaTeX for tick labels and set font size
end

% Save the figure
savefig('evol_modes_cel.fig');
