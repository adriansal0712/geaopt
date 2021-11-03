function plot_iter_curve(best_hist_gea, best_hist_pso, best_hist_bat, Num_iterations, func_name)

semilogy(best_hist_gea(:,1), best_hist_gea(:,2), 'r.-');
hold on
semilogy(best_hist_pso(:,1), best_hist_pso(:,2), 'b.-');
hold on
semilogy(best_hist_bat(:,1), best_hist_bat(:,2), 'g.-');
hold on
axis([0 Num_iterations -20 1000]);
switch func_name
    case 'F1' 
        title('Error results', 'FontSize', 16); 
        
    case 'F2'
        title('Error results', 'FontSize', 16); 
        
    case 'F3' 
        title('Error results', 'FontSize', 16); 
        
    case 'F4'
        title('Error results', 'FontSize', 16); 
        
    case 'F5' 
        title('Error results', 'FontSize', 16); 
        
    case 'F6'
        title('Error results', 'FontSize', 16);       
end
xlabel('no of iterations', 'FontSize', 16);
ylabel('minima', 'FontSize', 16);
legend('GEA optimization', 'PSO Optimization', 'BAT Algorithm', 'FontSize', 16);
hold off
end