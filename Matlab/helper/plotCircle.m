%% plot a black circle
function [ h ] = plotCircle( r, x,y )  
    display(['ploting circle: r:', num2str(r) , ', x:', num2str(x) , ', y:', num2str(y)]);

    th = 0:pi/50:2*pi;
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    h = plot(xunit, yunit,'color','black','LineWidth',2);

end

