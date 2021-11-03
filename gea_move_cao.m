% move all particles toward global minimum
% This method is based on Cao's method, (not GEA, but using frequency tuning)
function [best, fmin_gea, best_hist] = gea_move_notlive_cao(pars,options,coordinates)

Qmin=0;  %% lower limits of frequency
Qmax=2;  %% upper limits of frequency

% Obtain coordinates

if nargin<3
    coordinates.x = [];
    coordinates.y = [];
    coordinates.z = [];
end

xplot = coordinates.x;
yplot = coordinates.y;
zplot = coordinates.z;

figure

Num_iterations = options.Num_iterations;    % number of iterations
Tmax = Num_iterations;
n = options.n;                              % no of solutions
func_name = pars.fgname;                    % func_name
d = pars.nvar;                              % nvar of pars 

Q=zeros(n,1);                                   % frequency - for Cao's method only

[Sol, Fitness, fobj, lb,ub] = initpar(pars,options);   % Call init_gea function to 
                                                % initialize solutions

                                                
% Find current best
[fmin_gea, I] = min(Fitness);
best = Sol(I,:);


for iter = 1: Num_iterations
    
    switch func_name
        case 'F1' 
            
            hold off;
            contour(xplot,yplot,zplot, 15); 
            hold on;
            title('De Jong sphere function using GEA', 'FontSize', 16); 
            plot(Sol(:,1), Sol(:,2), 'b.', best(1), best(2),'r*');
            axis([lb(1) ub(1) lb(2) ub(2)]);


        case 'F2'
            hold off;
            contour(xplot,yplot,zplot, 15); 
            hold on;
            title('Schwefel 2.2 function using GEA', 'FontSize', 16); 
            plot(Sol(:,1), Sol(:,2), 'b.', best(1), best(2),'r*');
            axis([lb(1) ub(1) lb(2) ub(2)]);

        case 'F3' 
            hold off;
            contour(xplot,yplot,zplot, 15); 
            hold on;
            title('Griewangk function using GEA', 'FontSize', 16); 
            plot(Sol(:,1), Sol(:,2), 'b.', best(1), best(2),'r*');
            axis([lb(1) ub(1) lb(2) ub(2)]);

        case 'F4'
            hold off;
            contour(xplot,yplot,zplot, 15); 
            hold on;
            title('Rosenbrock function using GEA', 'FontSize', 16); 
            plot(Sol(:,1), Sol(:,2), 'b.', best(1), best(2),'r*');
            axis([lb(1) ub(1) lb(2) ub(2)]);

        case 'F5' 
            hold off;
            contour(xplot,yplot,zplot, 15); 
            hold on;
            title('Rastrigin function using GEA', 'FontSize', 16); 
            plot(Sol(:,1), Sol(:,2), 'b.', best(1), best(2),'r*');
            axis([lb(1) ub(1) lb(2) ub(2)]);

        case 'F6'
            hold off;
            contour(xplot,yplot,zplot, 15); 
            hold on;
            title('Michalewicz function using GEA', 'FontSize', 16); 
            plot(Sol(:,1), Sol(:,2), 'b.', best(1), best(2),'r*');
            axis([lb(1) ub(1) lb(2) ub(2)]);
            
        case 'F7'
            hold off;
            contour(xplot,yplot,zplot, 50); 
            hold on;
            title('Michalewicz function using GEA', 'FontSize', 16); 
            plot(Sol(:,1), Sol(:,2), 'b.', best(1), best(2),'r*');
            axis([lb(1) ub(1) lb(2) ub(2)]);
            
        case 'F8'
            hold off;
            contour(xplot,yplot,zplot, 50); 
            hold on;
            title('Non-collocated control using GEA', 'FontSize', 16); 
            plot(Sol(:,1), Sol(:,2), 'b.', best(1), best(2),'r*');
            axis([lb(1) ub(1) lb(2) ub(2)]);
            
        otherwise
            scatter(iter,fmin_gea,'b.');
            hold on
            title('Cost Function');
            xlim([0 Tmax]);
            xlabel("Number of iterations");
            ylabel("Cost Function");
            
    end
    
    
    
    for i = 1:n,

        S(i,:) = Sol(i,:);

        % 1 - Crossover operation - using Cao's method
        % Check condition: new position must be within bounds

        % Using Cao's method - using frequency tuning

        Q(i)=Qmin+(Qmax-Qmin)*rand;
        S(i,:)=Sol(i,:)+(best-Sol(i,:)).*Q(i);

        % Simple bounds
        Flag4up=S(i,:)>ub;
        Flag4low=S(i,:)<lb;
        S(i,:)=S(i,:).*(~(Flag4up+Flag4low))+ub.*Flag4up+lb.*Flag4low;

        % 2 - Mutation operator
        c = 0.2;                            % limiting parameter
        p = c*log(Tmax/(Tmax-iter));                    %   mutation probability
        M = 0.5*(ub - lb);         % Cao uses mutation vector = 20% of range

        % Mutation using total parameters
        if rand < p
            S(i,:)=S(i,:)+M.*unifrnd(-1,1,[1,d]); 
            Flag4up=S(i,:)>ub;
            Flag4low=S(i,:)<lb;
            S(i,:)=S(i,:).*(~(Flag4up+Flag4low))+ub.*Flag4up+lb.*Flag4low;
        end 


        % 3 - Local Search using Cao's method

        L = 0.1.*(ub - lb);        % Local search vector
        if rand < p
           S(i,:)=best+L.*unifrnd(-1,1,[1,d]);
            % Simple bounds
            Flag4up=S(i,:)>ub;
            Flag4low=S(i,:)<lb;
            S(i,:)=S(i,:).*(~(Flag4up+Flag4low))+ub.*Flag4up+lb.*Flag4low;

        end

        % Select new individual. If objective function is lower than in
        % previous iteration, select individual

        if fobj(S(i,:),pars) < Fitness(i),
                Fitness(i) = fobj(S(i,:),pars);
                Sol(i,:) = S(i,:);
        end
            
            % if minimum value is lower than the global minimum, select individual
        if Fitness(i) < fmin_gea,
            best = Sol(i,:);
            fmin_gea = Fitness(i);
        end

    end
    
    drawnow;
    
    best_hist(iter,:) = [iter, fmin_gea];

end

end