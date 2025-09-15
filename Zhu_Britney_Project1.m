% AMS 595 Project 1: How to find pi?

% Question 1: Using For Loops to Compute 

% Trying 1000 random points
N = 1e3
true_pi = pi %comparing estimated with true value of pi
x = rand(N,1) % random values of X from 0 to 1
y = rand(N,1) % random values of Y from 0 to 1

inside = (x.^2 + y.^2) <= 1;  % inside if less than or equal to 1
pi_est = 4 * sum(inside) / N; % formula to find calculation of pi

fprintf('Estimated pi with %d points = %.5f\n', N, pi_est); % prints the estimate of pi

%plots the estimation of pi of the 1000points compared to true value of pi
figure;
subplot(1,2,1);
plot(N, pi_est, '-o');
hold on;
yline(true_pi,'r--','True \pi');
xlabel('Number of Points');
ylabel('Estimated \pi');
title('The Estimation of \pi (For Loop)');

% Trying different fixed random points.
%% Initialize variables
N_lengths = [1e2, 1e3, 1e4, 1e5, 1e6];   % different point counts
pi_true = pi; %true estimate of pi
pi_estimates = zeros(size(N_lengths)); % estimate arrays to hold values of what we got but setting it to 0 value.
errors = zeros(size(N_lengths)); % errors arrays to hold values of what we got but setting it to 0 value.
times = zeros(size(N_lengths)); % times arrays to hold values of what we got but setting it to 0 value.

for k = 1:length(N_lengths) %for each of the fixed random time points
    N = N_lengths(k);
    tic;  % start the executation timer
    
    inside = 0;
    for i = 1:N
        x = rand(); % random values of X from 0 to 1
        y = rand(); % random values of Y from 0 to 1
        if x^2 + y^2 <= 1 % to check if inside circle, +1 inside if it is.
            inside = inside + 1;
        end
    end
    
    pi_est = 4 * inside / N; %formula to calc pi
    times(k) = toc;  % stop the execuation timer
    
    pi_estimates(k) = pi_est;
    errors(k) = abs(pi_true - pi_est); % see the diff in error
end

% Plotting the results, first graph to show the estimation of pi, second to
% show the error and time cost
figure;
subplot(1,2,1);
plot(N_lengths, pi_estimates, '-o');
hold on;
yline(pi_true,'r--','True \pi');
xlabel('Number of Points');
ylabel('Estimated \pi');
title('The Estimation of \pi (For Loop)');

subplot(1,2,2);
plot(times, errors, '-o');
xlabel('Execution Time (s)');
ylabel('Error |π - π_{est}|');
title('Error vs Time Cost');


% Question 2: Using While loops

% To compute pi to a certain number of significant figures using while
% loops we used 2, 3 and 4 as examples. Cannot use the true pi value.
% The gameplan is similar to for loops but to work backwards.
sig_level = [2, 3, 4];  
results = struct; % stores the results for each precision level

% Loop through each precision level
for p = 1:length(sig_level)
    
    digits_required = sig_level(p);   % how many sig figs we want
    
    % We can build tolerance to have a small margin of error that we can
    % accept between 2 values. This can help prevent endless loop and saves
    % more estimation time.
    tolerance = 0.5 * 10^(-digits_required); % we do tolerance base on current sig.level, so if it is 2 sig, then its 0.5*10^-2
    
    % Initializing the variables
    points_inside_circle = 0;
    total_points = 0;
    pi_estimate = 0;
    pi_previous = 0;
    
    % Now we use the while loops to make points until the estimate gets
    % close to the signifance level wanted
    while true
        % Again, we generate a random point (x,y) in 0 <= x,y <= 1
        x = rand();
        y = rand();
        
        % Count total points to see number of iterations
        total_points = total_points + 1;
        
        % USing te formula of circle to see if the point counts inside.
        if x^2 + y^2 <= 1
            points_inside_circle = points_inside_circle + 1;
        end
        
        % Compute current estimate of pi using same formula from Q1
        pi_previous = pi_estimate; % store the old estimate
        pi_estimate = 4 * points_inside_circle / total_points;
        
        % Stop if estimate has stop changing significantly (as our difference gets smaller than our tolerance)
        % and at least100 points to avoid very early stopping and more accurate calculation
        if total_points > 100 && abs(pi_estimate - pi_previous) < tolerance
            break;
        end
    end
    
    % Store the final results
    results(p).digits = digits_required;
    results(p).iterations = total_points;
    results(p).pi_est = pi_estimate;
    
    % Display the final results based on the different significant figures
    fprintf('%d significant figures: π ≈ %.*f after %d iterations\n', ...
            digits_required, digits_required, pi_estimate, total_points);
end

% Question 3: Creating a function with visual graphs for our pi calculation.


% Example: compute pi to 4 significant figures
pi_est_final = PiEstimation(4);


% Function Definition

function pi_estimate = PiEstimation(sig_figure)
    % This function uses Monte Carlo algorithm to estimate the calculation
    % of pi by using while loops with the given input of significance
    % figures that the user inputs
    % sig_figure: number of significant figures the user inputs

    % Creates a tolerance to prevent endless loops
    tolerance = 0.5 * 10^(-sig_figure);

    % Initializing the variables again
    points_inside_circle = 0;
    total_points = 0;
    pi_previous = 0;
    pi_estimate = 0;

    % Create figure for plotting points
    figure;
    hold on;
    axis equal;
    axis([0 1 0 1]); %creates a one by one figure
    xlabel('x'); ylabel('y');
    title(['Estimation of \pi to ', num2str(sig_figure), ' sig figs']);

    % Now we use the while loops to make points until the estimate gets
    % close to the signifance level wanted
    while true
        x = rand();
        y = rand();
        total_points = total_points + 1;

        % Check if point is inside circle
        if x^2 + y^2 <= 1
            points_inside_circle = points_inside_circle + 1;
            plot(x, y, 'g.'); % visually can see green for points inside
        else
            plot(x, y, 'r.'); % visually can see red for points outside
        end

        % Update the estimate calculated
        pi_previous = pi_estimate;
        pi_estimate = 4 * points_inside_circle / total_points;

        % Stop if estimate is less than tolenerance but ahave at least 100
        % points
        if total_points > 100 && abs(pi_estimate - pi_previous) < tolerance
            break;
        end

        % Update plot every 500 points (optional, can help with speeds up plotting as sig level increases)
        if mod(total_points,500) == 0
            drawnow;
        end
    end

    % Show final result of the pi calculated
    pi_str = sprintf(['pi ≈ %.', num2str(sig_figure), 'f'], pi_estimate);
    disp(['Final result: ', pi_str]);

    % Display and prints result on plot in blue font
    text(0.5, 0.5, pi_str, 'FontSize', 14, 'Color', 'b', ...
        'HorizontalAlignment','center');
end