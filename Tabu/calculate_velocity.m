%This function calculate the wind velocity deficit percentage resulted from
%the wake loss model def
%assuming wake spreading constant is 2
%radius of turbine is 20m

function vel = calculate_velocity(x, y) 
    global gridSize matrixSize matrix windVel
    
    gridSize = 80;
    matrixSize=100;
    %thrus coefficient of the turbine
    ct = 0.88;
    k = 2;
    R = 20;
    
    vel_def_total = 0;
    for i = 1 : x-1
        for j = 1 : matrixSize
            if matrix(i,j)==1
                if check_wake(i*gridSize, j*gridSize, x*gridSize, y*gridSize)==1
                    numerator
                    vel_def_cur = (1-sqrt(1-ct))/(1+k*(x-i)*gridSize/R)^2;
                    vel_def_total = vel_def_total+vel_def_cur^2;
                end
            end
        end
    end
    vel_def = sqrt(vel_def_total);
    vel = windVel * vel_def;
end

%check if turbine 2(x2,y2) is in wake of turbine 1(x1,y1)
function check = check_wake(x1,y1,x2,y2)
    %wind direction is 0 degree
    theta =0;
    alpha = arctan(2);
    R=20;
    k=2;
    
    numerator = (x2-x1)*cos(theta)+(y2-y1)*sin(theta)+R/k;
    denominator = sqrt( (x2-x1+R/k*cos(theta))^2 + (y2-y1+R/k*cos(theta))^2);
    beta = arccos( numerator/denominator );
    if(beta<alpha)
        check=1;
    else
        check=0;
    end
end