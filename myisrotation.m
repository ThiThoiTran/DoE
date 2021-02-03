% myisrotation of two given binary configurations
function [r] = myisrotation(y,y1)
n = length(y);
if n~= length(y1)
    fprintf('check again the dimension');
end
r = 0;
for i = 1:n
    if isequal(y,circshift(y1,[0,i]))
        r = 1;
        break;
    end
end
end