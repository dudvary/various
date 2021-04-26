%% Plot either split violin or entire violin
% Input: 
% - data: vector with data [1 x N]
% - xpos: Position on x axis where violin is split (default: 1)
% - side: 'r' or 'l' for side of split violin; 'b' for normal
%   violin (default: 'b'); 
% - FaceColor: Color of Area under Violin (default: 'k'); 
function violin_split(data,xpos,side,FaceColor)

    if nargin<4
        FaceColor = 'b'; 
    elseif nargin <3
        side = 'b'; 
    elseif nargin <2
        xpos = 1; 
    end

    EdgeColor = 'k'; 
    Alpha = 0.9; 
    ColorMedian = 'k';    
    wdth = 0.4; 
    
    [F, U, ~]=ksdensity(data);
    F=F/max(F)*wdth; %normalize

    MED = nanmedian(data); % Median
    Perc25 = prctile(data,25); % First Precentile
    Perc75 = prctile(data,75); % First Precentile
    
    if strcmp(side,'l');
        h=fill([xpos*ones(size(F));flipud(xpos-F)]',[U;flipud(U)]',FaceColor,'FaceAlpha',Alpha,'EdgeColor',EdgeColor,'LineWidth',1);
        hold on
        plot([interp1(U,xpos*ones(size(F)),MED), interp1(flipud(U),flipud(xpos-F),MED) ],[MED MED],'-','Color',ColorMedian,'LineWidth',2);
        hold on
        plot([interp1(U,xpos*ones(size(F)),Perc25), interp1(flipud(U),flipud(xpos-F),Perc25) ],[Perc25 Perc25],'-','Color',ColorMedian,'LineWidth',1);
        hold on
        plot([interp1(U,xpos*ones(size(F)),Perc75), interp1(flipud(U),flipud(xpos-F),Perc75) ],[Perc75 Perc75],'-','Color',ColorMedian,'LineWidth',1);
        hold on;
    elseif strcmp(side,'r');
        h=fill([F+xpos;flipud(xpos*ones(size(F)))]',[U;flipud(U)]',FaceColor,'FaceAlpha',Alpha,'EdgeColor',EdgeColor,'LineWidth',1);
        hold on
        plot([interp1(U,F+xpos,MED), interp1(flipud(U),flipud(xpos*ones(size(F))),MED) ],[MED MED],'-','Color',ColorMedian,'LineWidth',2);
        hold on
        plot([interp1(U,F+xpos,Perc25), interp1(flipud(U),flipud(xpos*ones(size(F))),Perc25) ],[Perc25 Perc25],'-','Color',ColorMedian,'LineWidth',1);
        hold on
        plot([interp1(U,F+xpos,Perc75), interp1(flipud(U),flipud(xpos*ones(size(F))),Perc75) ],[Perc75 Perc75],'-','Color',ColorMedian,'LineWidth',1);
        hold on; 
    elseif strcmp(side,'b')
        h=fill([F+xpos fliplr(xpos-F)]',[U fliplr(U)]',FaceColor,'FaceAlpha',Alpha,'EdgeColor',EdgeColor);
        hold on
        plot([interp1(U,F+xpos,MED), interp1(flipud(U),flipud(xpos-F),MED) ],[MED MED],'-','Color',ColorMedian,'LineWidth',2);
        hold on
        plot([interp1(U,F+xpos,Perc25), interp1(flipud(U),flipud(xpos-F),Perc25) ],[Perc25 Perc25],'-','Color',ColorMedian,'LineWidth',1);
        hold on
        plot([interp1(U,F+xpos,Perc75), interp1(flipud(U),flipud(xpos-F),Perc75) ],[Perc75 Perc75],'-','Color',ColorMedian,'LineWidth',1);
    else
       warning(['side needs to be either r or l']);  
       error([side ' is not a valid input!']);
    end
end