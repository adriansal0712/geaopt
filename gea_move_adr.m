%% GEA Move: function for moving particles from the initial location to the new location
% Reference: A Guiding Evolutionary Algorithm for Optimization
% move all particles toward global minimum

function [x_gea, f_gea, best_hist] = gea_move_adr(pars,options,coordinates)

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
[f_gea, I] = min(Fitness);
x_gea = Sol(I,:);

for iter = 1: Num_iterations
    
    
    switch func_name
        case 'F1' 
            
            hold off;
            contour(xplot,yplot,zplot, 15); 
            hold on;
            title('De Jong sphere function using GEA', 'FontSize', 16); 
            plot(Sol(:,1), Sol(:,2), 'b.', x_gea(1), x_gea(2),'r*');
            axis([lb(1) ub(1) lb(2) ub(2)]);


        case 'F2'
            hold off;
            contour(xplot,yplot,zplot, 15); 
            hold on;
            title('Schwefel 2.2 function using GEA', 'FontSize', 16); 
            plot(Sol(:,1), Sol(:,2), 'b.', x_gea(1), x_gea(2),'r*');
            axis([lb(1) ub(1) lb(2) ub(2)]);

        case 'F3' 
            hold off;
            contour(xplot,yplot,zplot, 15); 
            hold on;
            title('Griewangk function using GEA', 'FontSize', 16); 
            plot(Sol(:,1), Sol(:,2), 'b.', x_gea(1), x_gea(2),'r*');
            axis([lb(1) ub(1) lb(2) ub(2)]);

        case 'F4'
            hold off;
            contour(xplot,yplot,zplot, 15); 
            hold on;
            title('Rosenbrock function using GEA', 'FontSize', 16); 
            plot(Sol(:,1), Sol(:,2), 'b.', x_gea(1), x_gea(2),'r*');
            axis([lb(1) ub(1) lb(2) ub(2)]);

        case 'F5' 
            hold off;
            contour(xplot,yplot,zplot, 15); 
            hold on;
            title('Rastrigin function using GEA', 'FontSize', 16); 
            plot(Sol(:,1), Sol(:,2), 'b.', x_gea(1), x_gea(2),'r*');
            axis([lb(1) ub(1) lb(2) ub(2)]);

        case 'F6'
            hold off;
            contour(xplot,yplot,zplot, 15); 
            hold on;
            title('Michalewicz function using GEA', 'FontSize', 16); 
            plot(Sol(:,1), Sol(:,2), 'b.', x_gea(1), x_gea(2),'r*');
            axis([lb(1) ub(1) lb(2) ub(2)]);
            
        case 'F7'
            hold off;
            contour(xplot,yplot,zplot, 50); 
            hold on;
            title('Non-collocated control using GEA', 'FontSize', 16); 
            plot(Sol(:,1), Sol(:,2), 'b.', x_gea(1), x_gea(2),'r*');
            axis([lb(1) ub(1) lb(2) ub(2)]);
            
        case 'F8'
            hold off;
            contour(xplot,yplot,zplot, 50); 
            hold on;
            title('Non-collocated control using GEA', 'FontSize', 16); 
            plot(Sol(:,1), Sol(:,2), 'b.', x_gea(1), x_gea(2),'r*');
            axis([lb(1) ub(1) lb(2) ub(2)]);
            
        otherwise
            scatter(iter,f_gea,'b.');
            hold on
            title('Cost Function');
            xlim([0 Tmax]);
            xlabel("Number of iterations");
            ylabel("Cost Function");
            
    end
    
%     % Using for loop 
%     for i = 1:n,
%         % 1 - Crossover operation - using original GEA method
%         % Check condition: new position must be within bounds
%         S(i,:) = Sol(i,:);
%         beta = 2*rand;
%         S(i,:)=Sol(i,:)+(x_gea-Sol(i,:))*beta;
%         
%         
%         % Simple bounds
%         Flag4up=S(i,:)>ub;
%         Flag4low=S(i,:)<lb;
%         S(i,:)=S(i,:).*(~(Flag4up+Flag4low))+ub.*Flag4up+lb.*Flag4low;
%         
%         % 2 - Mutation operator
%         % Mutation using individual parameters
%         c = 0.2;                            % limiting parameter
%         p = c*log(Tmax/(Tmax-iter));                    %   mutation probability
% 
%         Mx = max(S(i,:)-lb, ub - S(i,:));
%         
%         r = rand(1,d)<p;
%         
%         S(i,:) = S(i,:).*(~r) + (S(i,:) + unifrnd(-1,1,[1 d]).*Mx).*(r); % applying mutation to individual dimensions only
%        
%         % Simple bounds
%         Flag4up=S(i,:)>ub;
%         Flag4low=S(i,:)<lb;
%         S(i,:)=S(i,:).*(~(Flag4up+Flag4low))+ub.*Flag4up+lb.*Flag4low;
%         
%         % 3 - Local Search
%         % Local search using individual parameters
%         L = 0.1.*(ub - lb);        % Local search vector
% 
%         % Local search using individual parameters
%         r = rand(1,d)<p;
% 
%         S(i,:) = S(i,:).*(~r) + r.*(S(i,:) + unifrnd(-1,1,[1 d]).*L); % applying local search to individual dimensions only
% 
%         % Simple bounds
%         Flag4up=S(i,:)>ub;
%         Flag4low=S(i,:)<lb;
%         S(i,:)=S(i,:).*(~(Flag4up+Flag4low))+ub.*Flag4up+lb.*Flag4low;
%         
%         
%         % Selection
%         if fobj(S(i,:),pars) < Fitness(i),
%             Fitness(i) = fobj(S(i,:),pars);
%             Sol(i,:) = S(i,:);
%         end
%             
%         % if minimum value is lower than the global minimum, select individual
%         if Fitness(i) < f_gea,
%             x_gea = Sol(i,:);
%             f_gea = Fitness(i);
%         end
%         
%     end
%     
%     drawnow;
%     
%     best_hist(iter,:) = [iter, f_gea];
%     x_gea
%     f_gea
% end
%     
        
% Vectorized method        
    S = Sol;

% % 1 - Crossover operation - using original GEA method
% % Check condition: new position must be within bounds
    
    beta = 2*rand([n,1]);                                          % step length of position increment
    S = Sol + (x_gea-Sol).*beta;

% %     Simplebounds
    Flag4up=S>ub;
    Flag4low=S<lb;
    S=S.*(~(Flag4up+Flag4low))+ub.*Flag4up+lb.*Flag4low;

% % 2 - Mutation operator
    c = 0.2;                            % limiting parameter
    p = c*log(Tmax/(Tmax-iter));                    %   mutation probability

    Mx = max(S-lb, ub - S);

% %     Mutation using individual parameters             
    r = rand(n,d)<p;
    S = S.*(~r) + (S + unifrnd(-1,1,[1 d]).*Mx).*(r);
    Flag4up=S>ub;
    Flag4low=S<lb;
    S=S.*(~(Flag4up+Flag4low))+ub.*Flag4up+lb.*Flag4low;

% % 3 - Local Search
% %   Local search using individual paramameters
        
    r = rand(n,d)<p;
    L = 0.1.*(ub - lb);        % Local search vector
    r = rand(1,d)<p;
    S = S.*(~r) + (S + unifrnd(-1,1,[n d]).*L).*(r);
               
    Flag4up=S>ub;
    Flag4low=S<lb;
    S=S.*(~(Flag4up+Flag4low))+ub.*Flag4up+lb.*Flag4low;

% % 4 - Selection
    for i = 1:n
%             Fitness(i) = fobj(S(i,:),pars);

            fitnessnew = fobj(S(i,:),pars);
            if fitnessnew < Fitness(i),
                Fitness(i) = fitnessnew;
                Sol(i,:) = S(i,:);
            end
            
            % if minimum value is lower than the global minimum, select individual
            if Fitness(i) < f_gea,
                x_gea = Sol(i,:);
                f_gea = Fitness(i);
            end
            
    end

    drawnow;
    
    best_hist(iter,:) = [iter, f_gea];
%     hold off;

    
end

% END of VECTORIZED METHOD
          

%     end 
            

    % End of Mutation using individual parameters

    % Mutation using total parameters - Comment out if not in use
%     if rand < p
%         S(i,[1 2])=S(i,[1 2])+Mx.*unifrnd(-1,1,[1,2]);
%         
%     end 
%     Flag4up=S(i,[1:2])>ub;
%     Flag4low=S(i,[1:2])<lb;
%     S(i,[1:2])=S(i,[1:2]).*(~(Flag4up+Flag4low))+ub.*Flag4up+lb.*Flag4low;

        
        % Applying local search to complete parameters
    %     if rand<p,
    %         S(i,[1 2]) = best([1 2]) + L.*unifrnd(-1,1,[1 2]);
    %     
    %         % Applying simple bounds
    % 
    %         Flag4up=S(i,[1:2])>ub;
    %         Flag4low=S(i,[1:2])<lb;
    %         S(i,[1:2])=S(i,[1:2]).*(~(Flag4up+Flag4low))+ub.*Flag4up+lb.*Flag4low;
    %     end

        % End of Local search


        

%         for i = 1:n
%             Fitness(i) = fobj(S(i,:),pars);
            
            
%             if fobj(S(i,:),pars) < Fitness(i),
%                 Fitness(i) = fobj(S(i,:),pars);
%                 Sol(i,:) = S(i,:);
%             end
%             
%             % if minimum value is lower than the global minimum, select individual
%             if Fitness(i) < f_gea,
%                 x_gea = Sol(i,:);
%                 f_gea = Fitness(i);
%             end
            
%         end
          
%         % Find current best
%         [f_gea, I] = min(Fitness);
%         x_gea = Sol(I,:);
        

%     end 
    


end