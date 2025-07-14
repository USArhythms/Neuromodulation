function f_plotBar(varargin)

data = varargin{1};

p = inputParser;
addParameter(p,'colors',[0, 0, 0]);
addParameter(p,'ylabel','');
addParameter(p,'legend','');
addParameter(p,'title','');

parse(p,varargin{2:end});

plotLegend = p.Results.legend;
colors = p.Results.colors;

if iscell(plotLegend)
    plotLegend = [plotLegend; repmat({''},1,numel(plotLegend))];
    plotLegend = plotLegend(:);
end

N = numel(data);
dataMean = zeros(N,1);
dataSEM = zeros(N,1);
for i = 1:N
    dataMean(i) = mean(data{i});
    dataSEM(i) = std(data{i},0)/sqrt(numel(data{i}));
end

hold on;
cIdx = round(linspace(1,size(colors,1),numel(data)));

for i = 1:N
    b(i) = bar(i,dataMean(i));
    scatter(i*ones(numel(data{i}),1),data{i},70,'filled',XJitter='randn',XJitterWidth=0.3,MarkerFaceColor=colors(cIdx(i),:),MarkerFaceAlpha=0.5);
    b(i).FaceColor = 'flat';
    b(i).CData = [1 1 1];
    b(i).ShowBaseLine = 'off';
    b(i).EdgeColor = colors(cIdx(i),:);
    b(i).LineWidth = 3;
end

er = errorbar(1:N,dataMean,dataSEM);

er.Color = [0 0 0];
er.LineStyle = 'none';
er.LineWidth = 3;

ax = gca;
ax.XAxis.Visible = 'off';

ylabel(p.Results.ylabel);
legend(plotLegend);
title(p.Results.title);
set(gca,'FontSize',14);

end