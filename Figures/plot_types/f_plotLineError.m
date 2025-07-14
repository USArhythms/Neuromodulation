function f_plotLineError(x,y,error,color)

x = x(:);
y = y(:);
error = error(:);

if nargin < 4
    color = get(groot,'defaultAxesColorOrder');
    color = color(1,:);
end

hold on;
fill([x;flipud(x)],[(y+error);flipud(y-error)],color,FaceAlpha=0.3,EdgeColor='none');
plot(x,y,Color=color,LineWidth=2);

end