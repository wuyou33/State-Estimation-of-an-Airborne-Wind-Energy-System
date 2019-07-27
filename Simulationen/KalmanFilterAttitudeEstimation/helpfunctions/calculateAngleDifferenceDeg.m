function diff = calculateAngleDifferenceDeg(angle1, angle2)
    % Problem: if one angle becomes greater than 180° ther is a
    % wrap around to -180° (pi rad). In this short time the angle
    % difference is 360°!              
    diff = 180 - abs(abs(angle1 - angle2) - 180);    
    diff(angle1<angle2) = diff(angle1<angle2) * -1;
end