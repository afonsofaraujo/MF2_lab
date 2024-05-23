function integral_approx = int_trap(x, y, a, b, plot_flag)
    idx_a = find(x >= a, 1, 'first');  % Index of first x value >= a
    idx_b = find(x <= b, 1, 'last');   % Index of last x value <= b
    x_integrate = x(idx_a:idx_b);
    y_integrate = y(idx_a:idx_b);
    h = diff(x_integrate);  % Differences between consecutive x values
    integral_approx = sum(h .* (y_integrate(1:end-1) + y_integrate(2:end)) / 2);
    if plot_flag
        figure;
        plot(x, y, 'b.-');  % Plot the data points
        hold on;
        plot(x_integrate, y_integrate, 'r.-');  % Highlight the interval for integration
        for k = 1:length(x_integrate) - 1
            fill([x_integrate(k), x_integrate(k+1), x_integrate(k+1), x_integrate(k)], ...
                 [0, 0, y_integrate(k+1), y_integrate(k)], 'r', 'FaceAlpha', 0.3);
        end
        title('Trapezoid Rule for Numerical Integration');
        xlabel('x');
        ylabel('y');
        legend('Data Points', 'Interval for Integration', 'Trapezoids', 'Location', 'northwest');
        set(gcf, 'Color', 'w');
        grid on;
        hold off;
    end
end
