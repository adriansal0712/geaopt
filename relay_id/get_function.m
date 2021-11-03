function [G, relay, Tds, num, den, Tsim, tol] = get_function(func_name, s)

switch func_name
    case 'P1'   % P1 = 1/((s+1)*(0.1*s+1)*(0.01*s+1)*(0.001*s+1))
        num = 1;
        den = sym2poly((s+1)*(0.1*s+1)*(0.01*s+1)); % *(0.001*s+1)
        assignin('caller','G',tf(num, den));
        assignin('caller','Tds',0);
        assignin('caller','relay',[0.1 -0.1 2 -1]);
        assignin('caller','Tsim',20);
        assignin('caller','tol',0.05);
        
    case 'P2'   % P2 = 1/(s+1)^4
        
        assignin('base','num',1);
        assignin('base','den',sym2poly((s+1)^4)); % *(0.001*s+1)
        assignin('base','G',tf(num, den));
        assignin('base','Tds',0);
        assignin('base','relay',[0.1 -0.1 2 -1]);
        assignin('base','Tsim',20);
        assignin('base','tol',0.1);
        
    case 'P3'    % P3 = e^-s/(0.05*s+1)^2
        
        assignin('base','num',1);
        assignin('base','den',sym2poly((0.05*s+1)^2)); % *(0.001*s+1)
        assignin('base','G',tf(num,den, 'InputDelay', 1));
        assignin('base','Tds',0);
        assignin('base','relay',[0.1 -0.1 2 -1]);
        assignin('base','Tsim',20);
        assignin('base','tol',0.5);
        
    case 'P4' % P4 = 1/0.5*s^2 + s + 1)^2
        
        assignin('base','num',1);
        assignin('base','den',sym2poly((0.5*s^2 + s + 1)^2)); % *(0.001*s+1)
        assignin('base','G',tf(num,den));
        assignin('base','Tds',0);
        assignin('base','relay',[0.1 -0.1 2 -1]);
        assignin('base','Tsim',20);
        assignin('base','tol',0.5);
        
    case 'P5'   % P5 = (e^-s)/((s+1)*(0.3s+1)^2))
        
        assignin('base','num',1);
        assignin('base','den',sym2poly((s+1)*(0.3*s+1)^2)); % *(0.001*s+1)
        assignin('base','Tds',1);
        assignin('base','G',tf(num,den),'InputDelay',Tds);
        assignin('base','relay',[0.1 -0.1 2 -1]);
        assignin('base','Tsim',30);
        assignin('base','tol',0.1);
        
        
    case 'P6'   % P6 = (e^-4s)/((s+1)*(0.3*s+1)^2)
        
        assignin('base','num',1);
        assignin('base','den',sym2poly((s+1)*(0.3*s+1)^2)); % *(0.001*s+1)
        assignin('base','Tds',4);
        assignin('base','G',tf(num,den),'InputDelay',Tds);
        assignin('base','relay',[0.1 -0.1 2 -1]);
        assignin('base','Tsim',30);
        assignin('base','tol',0.3);

        
    case 'P7'   % function F6 acc to Thesis (e^-s)/((5*s+1)^5) 
        
        assignin('base','num',1);
        assignin('base','den',sym2poly((5*s+1)^5)); % *(0.001*s+1)
        assignin('base','Tds',1);
        assignin('base','G',tf(num,den),'InputDelay',Tds);
        assignin('base','relay',[0.1 -0.1 2 -1]);
        assignin('base','Tsim',200);
        assignin('base','tol',0.3);
        
        
end
