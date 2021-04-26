clear all
close all
clc

%%
filename = ""; % path/to/ratings.csv
opts = delimitedTextImportOptions("NumVariables", 13);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Var1", "YourRating", "DateRated", "Title", "URL", ...
                        "TitleType", "IMDbRating", "Runtimemins", "Year", ...
                        "Genres", "NumVotes", "Var12", "Var13"];
opts.SelectedVariableNames = ["YourRating", "DateRated", "Title", "URL", ...
                        "TitleType", "IMDbRating", "Runtimemins", "Year", ...
                        "Genres", "NumVotes"];
opts.VariableTypes = ["string", "double", "datetime", "string", "double", ...
                        "categorical", "double", "double", "double", ...
                        "categorical", "double", "string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["Var1", "Title", "Var12", "Var13"], ...
                                            "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var1", "Title", "TitleType", "Genres", ...
                                "Var12", "Var13"], "EmptyFieldRule", "auto");
opts = setvaropts(opts, "DateRated", "InputFormat", "yyyy-MM-dd");
opts = setvaropts(opts, "URL", "TrimNonNumeric", true);
opts = setvaropts(opts, "URL", "ThousandsSeparator", ",");

% Import the data
tbl = readtable(filename, opts);

% Convert to output type
MyRating = tbl.YourRating;
DateRated = tbl.DateRated;
% Title = tbl.Title;
% URL = tbl.URL;
TitleType = tbl.TitleType;
IMDbRating = tbl.IMDbRating;
% Runtimemins = tbl.Runtimemins;
Year = tbl.Year;
Genres = tbl.Genres;
% NumVotes = tbl.NumVotes;

% Clear temporary variables
clear opts tbl

yearRated = year(DateRated);
[C,ia,ic] = unique(yearRated);
yearCounts = accumarray(ic,1); 

%% Movies vs. TV Shows
f1 = figure(1);
clf;

% Histogram Ratings
subplot(1,3,1);
h = histogram(MyRating,[0:10]+0.5, ...
            'FaceColor',[0.7 0.7 0.7],'EdgeColor','w','FaceAlpha',1);
hold on;
errorbar(mean(MyRating),max(h.Values(:))*1.05,std(MyRating),'horizontal', ...
                        'ko','MarkerSize',10,'MarkerFaceColor','k'); 
set(gca,'Box','off','TickDir','out','XLim',[0.5 10.5],'XTick',[1:10], ...
    'YTick',[0:100:200],'YMinorTick','on','YLim',[0 max(h.Values(:))*1.1]);
xlabel('my rating');
ylabel('frequency');
axis square;

% My rating vs. IMDb rating
subplot(1,3,2);
scatter(MyRating,IMDbRating,30,'MarkerFaceColor','k', ...
                        'MarkerEdgeColor','none','MarkerFaceAlpha',0.3);
hold on;
plot([1 10],[1 10],':','Color',[0.5 0.5 0.5]);
p = polyfit(MyRating,IMDbRating,1);
yfit = polyval(p,[min(MyRating) max(MyRating)]); 
plot([min(MyRating) max(MyRating)],yfit,'-','Color',[0.5 0.5 1],'LineWidth',1.5);

r = corr(MyRating,IMDbRating);
text(1.5,9.5,['R = ' num2str(r,'%.2f')]);
xlabel('my rating');
ylabel('avg. IMDb rating');
axis([0.7 10.3 0.7 10.3]);
set(gca,'Box','off','TickDir','out','YTick',[1:10],'XTick',[1:10]);
axis square;

% Time vs. MyRating
labelStr = cell(1,length(yearCounts)); 
for i = 1:length(yearCounts)
   labelStr{i} = sprintf('%d (n=%d)',C(i),yearCounts(i));  
end

subplot(1,3,3);
boxplot(MyRating,yearRated,'Symbol','+','Color',[0.7 0.7 0.7], ... 
        'PlotStyle','compact','Labels',labelStr);
hold on;
set(gca,'Box','off','TickDir','out');
ylim([0.7 10.3]);
xlabel('year rated');
ylabel('my rating');
axis square;