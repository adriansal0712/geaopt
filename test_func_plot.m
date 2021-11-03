% Source: https://www.mathworks.com/matlabcentral/fileexchange/68981-bat-optimization-algorithm
% author: Abhishek Gupta
% modified by: Adrian Saldanha


function [coordinates,range,func_min] = test_func_plot(func_name)
[range,dim,fobj,func_min]=test_function(func_name);

subplot(1,2,1);

Ngrid = 100;    % number of points on grid, dx dy = spacing between pts
dx = (range(2) - range(1))/Ngrid;
dy = (range(4) - range(3))/Ngrid;
xgrid = range(1):dx:range(2);
ygrid = range(3):dy:range(4);
[coordinates.x coordinates.y] = meshgrid(xgrid, ygrid);

switch func_name 
    case 'F1'          % De Jong's sphere function
        coordinates.z=fobj(coordinates.x,coordinates.y);
        surfc(coordinates.x,coordinates.y,coordinates.z);
        title('De Jong''s Sphere function', 'FontSize', 16);
        func_min = 0;
        
    case 'F2'           % Schwefel 2.2 function
        coordinates.z=fobj(coordinates.x,coordinates.y);
        surfc(coordinates.x,coordinates.y,coordinates.z);
        title('Schwefel 2.2 function', 'FontSize', 16);

    case 'F3'           % Griewangk's function
        coordinates.z=fobj(coordinates.x,coordinates.y);
        surfc(coordinates.x,coordinates.y,coordinates.z);
        title('Griewangk function', 'FontSize', 16);
        
    case 'F4'           % Rosenbrock's function
        coordinates.z=fobj(coordinates.x,coordinates.y);
        surfc(coordinates.x,coordinates.y,coordinates.z);
        title('Rosenbrocks function', 'FontSize', 16);
        
    case 'F5'           % Rastrigin's function
        coordinates.z=fobj(coordinates.x,coordinates.y);
        surfc(coordinates.x,coordinates.y,coordinates.z);
        title('Rastrigins function', 'FontSize', 16);
    
    case 'F6'           % Michalewicz function
        coordinates.z=fobj(coordinates.x,coordinates.y);
        surfc(coordinates.x,coordinates.y,coordinates.z);
        title('Double-loop, two dof optimization', 'FontSize', 16);
        
    case 'F7'
        coordinates.z=fobj(coordinates.x,coordinates.y);
        surfc(coordinates.x,coordinates.y,coordinates.z);
        title('Non-collocated control optimization','FontSize',16);
                      
end    

subplot(1,2,2);
contour(coordinates.x,coordinates.y,coordinates.z, 15);


end