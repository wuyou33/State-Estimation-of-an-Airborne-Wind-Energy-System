% % https://de.mathworks.com/help/matlab/ref/uislider.html    

%function SliderFigure(time, p, gamma_inertia)  
global achsXx;
global achsXy;
global achsXz;
global achsYx;
global achsYy;
global achsYz;
global achsZx;
global achsZy;
global achsZz;

coordSystem = CoordSystem(1);
pos = p;
orientation = gamma_inertia;

fWidth = 500;
fHeight = 500;
fig = figure('Name', 'Position and Orientation','Position',[100 100 fWidth fHeight]);
figure(fig);
xlabel('X');
ylabel('Y');
zlabel('Z');
hold on;
minX = min(pos(1,:)); maxX = max(pos(1,:));
if minX < maxX
    xlim([min(pos(1,:)), max(pos(1,:))]);
else
    if minX == 0
        xlim([-1, 1]);
    else
        xlim([minX-abs(minX)*0.1, maxX+abs(maxX)*0.1]);
    end
end
minY = min(pos(2,:)); maxY = max(pos(2,:));
if minY < maxY
    ylim([min(pos(2,:)), max(pos(2,:))]);
else
    if minY == 0
        ylim([-1, 1]);
    else
        ylim([minY-abs(minY)*0.1, maxY+abs(maxY)*0.1]);
    end
end
minZ = min(pos(3,:)); maxZ = max(pos(3,:));
if minZ < maxZ
    zlim([min(pos(3,:)), max(pos(3,:))]);
else
    if minZ == 0
        zlim([-1, 1]);
    else
       zlim([minZ-abs(minZ)*0.1, maxZ+abs(maxZ)*0.1]);  
    end

end
%draw complete position path
plot3(pos(1,:), pos(2,:), pos(3,:), 'DisplayName', 'Real position');


% initial orientation
gamma = [orientation(1,1);
         orientation(2,1);
         orientation(3,1)];
actualPos = pos(:, 1);
[achsX, achsY, achsZ] = coordSystem.editDeg(actualPos, gamma);
achsXx = achsX(1,:);
achsXy = achsX(2,:);
achsXz = achsX(3,:);
achsYx = achsY(1,:);
achsYy = achsY(2,:);
achsYz = achsY(3,:);
achsZx = achsZ(1,:);
achsZy = achsZ(2,:);
achsZz = achsZ(3,:);
p3 = plot3(achsX(1,:), achsX(2,:), achsX(3,:), 'r', 'DisplayName', 'X');
p3.XDataSource = 'achsXx';
p3.YDataSource = 'achsXy';
p3.ZDataSource = 'achsXz';
p3 = plot3(achsY(1,:), achsY(2,:), achsY(3,:), 'g', 'DisplayName', 'Y');
p3.XDataSource = 'achsYx';
p3.YDataSource = 'achsYy';
p3.ZDataSource = 'achsYz';
p3 = plot3(achsZ(1,:), achsZ(2,:), achsZ(3,:), 'b', 'DisplayName', 'Z');
p3.XDataSource = 'achsZx';
p3.YDataSource = 'achsZy';
p3.ZDataSource = 'achsZz';
ax = gca;
ax.YDir = 'reverse';
%linkdata on;
hold off;

 % slider is on a uifigure because cannot placed on a figure
figSlider = uifigure('Name', 'Slider','Position',[100 100 fWidth fHeight]);

height = 3;
width = fWidth-40;
x = 20;
y = 20;
sld = uislider(figSlider,'Position',[x, fHeight-height-y, width, height], ...
                    'ValueChangedFcn',@(sld,event) updateValues(sld,time, pos, orientation, coordSystem));
sld.Limits = [min(time), max(time)];
sld.Value = min(time);


% Create ValueChangedFcn callback
function updateValues(sld, time, pos, orientation, coordSystem)
    global achsXx;
    global achsXy;
    global achsXz;
    global achsYx;
    global achsYy;
    global achsYz;
    global achsZx;
    global achsZy;
    global achsZz;
    actual_time = sld.Value;
    index = find(time >= actual_time,1);
    gamma = orientation(:,index);
    actualPos = pos(:,index);

    [achsX, achsY, achsZ] = coordSystem.editDeg(actualPos, gamma);
    achsXx = achsX(1,:);
    achsXy = achsX(2,:);
    achsXz = achsX(3,:);
    achsYx = achsY(1,:);
    achsYy = achsY(2,:);
    achsYz = achsY(3,:);
    achsZx = achsZ(1,:);
    achsZy = achsZ(2,:);
    achsZz = achsZ(3,:);
    refreshdata;
    drawnow;
end
