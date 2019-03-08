function [cRet,cRet2,cRet3,cRet4] = WithAPCA(pMat,nReps,ns,name)

    % function [cRet] = segmentM(pMat,nFourier,nReps)
    % ---------------------------------------
    % pMat      = 1D matrix of values
    % 
    % nReps     = max segments to look for
    %             DEFAULT = 30
    %             (bigger number takes more time)
    % ns        = smooth by ns median for medfilt1
    %

    %% set default arguments
    

    %% Flip matrix to the correct orientation if necessary
    [m n] = size(pMat);
    if m>n
        pMat = pMat';
    end

    kn=pMat;
    mMat=medfilt1(kn,ns);
    fig1=figure('Visible','off');
    hold all;
    f='chr\d+|chrX';
    newname=regexp(name,f,'match');
    
    a=[1:length(mMat)];
    a=a';
    plot(a*0.1,pMat,'b','LineWidth',1); 
    plot(a*0.1,mMat,'m','LineWidth',2); 
    legend([{'Original data'} {'medfilt data'}]);
    set(gca,'fontsize',16);
    title(newname);
    xlabel('Genome size (Mb)');
    
    ylabel('log2(normalized FPKM)');
    fname=sprintf('%s_data_smooth_%d.eps',name,ns);
    print(fig1,fname,'-depsc');
    hold off;
    close(fig1);
    %% Run APCA for each value of n
    % APCA will attempt to find the best trnasformation of the data into 
    % linear segments
    nRange          = 2:nReps; 
    distPlotData    = zeros(length(nRange),2);
    nCnt            = 0;
    
    fig13=figure('Visible','off');
    hold all;
    
    for n = nRange
        % Run APCA and return X,Y coordinates of segments
        [dAPCA,xAPCA,yAPCA] = aPCAMe(mMat,n);

        % plot segments for this n
        

        % calculate the error of this fit using Euclidean distance
        nDist = dist([mMat;yAPCA]');
        plot(xAPCA*0.1,yAPCA,'-');
        cLeg = [{sprintf('n = %i; d = %4.2f\n',n,nDist(1,2))}];
        % add legend with n and error
        

        % use gradient descent to calculate the inflexion point of the error
        % curve. Save everything to a matrix to allow us to plot this at the end
        % NOTE: THIS CAN CERTAINLY BE IMPROVED.
        nCnt = nCnt+1;
        distPlotData(nCnt,:) = [n nDist(1,2)];
        
        fprintf('n = %i; err = %4.2f\n',n,nDist(1,2));
        
    end
    legend(cLeg);
    set(gca,'fontsize',20);
    fname=sprintf('%s_apca_all_segments_%d.eps',name,ns);
    print(fig13,fname,'-depsc');
    hold off;
    close(fig13);
    
    % distPlotData
    nPoints=length(nRange);
    firstPoint = distPlotData(1,:);
    linevec=distPlotData(end,:)-firstPoint;
    lineVecN = linevec / sqrt(sum(linevec.^2));
    
    vecFromFirst = bsxfun(@minus, distPlotData, firstPoint);
    length(vecFromFirst);

    scalarProduct = dot(vecFromFirst, repmat(lineVecN,nPoints,1), 2);
    length(scalarProduct);
    vecFromFirstParallel = scalarProduct * lineVecN;
    vecToLine = vecFromFirst - vecFromFirstParallel;
    length(vecToLine);
    distToLine = sqrt(sum(vecToLine.^2,2));
    
    % remove the segment less then 3
    distPlotData(distPlotData(:,1)<3) = 0
    distPlotData( ~ all(distPlotData,2), : ) = [];
    distPlotData
    
%# now all you need is to find the maximum
    [maxDist,idxOfBestPoint] = sort(distToLine,'descend');
    distToLine;
    nBest=distPlotData(idxOfBestPoint(1),1)-1;
    nBest2=distPlotData(idxOfBestPoint(2),1)-1;
    nBest3=distPlotData(idxOfBestPoint(3),1)-1;
    nBest4=distPlotData(idxOfBestPoint(4),1)-1;
   
%# fplot
    fig2=figure('Visible','off');
    hold all;
    plot(distPlotData(:,1),distPlotData(:,2),'o-');
    
    plot(distPlotData(idxOfBestPoint(1),1), distPlotData(idxOfBestPoint(1),2), 'or','LineWidth',3);
    legend([ {'all APCA'} {sprintf('Best fit %i',nBest)}]);
    set(gca,'fontsize',20);
    xlabel('# segments (n)');
    ylabel('error');
    title(newname);
    fname=sprintf('%s_error_first.eps',name);
    print(fig2,fname,'-depsc');
    hold off;
    close(fig2);
    %[M,I]=min(distPlotData(:,2));
    
    nBest
    
    fig5=figure('Visible','off');
    hold all;
    plot(distPlotData(:,1),distPlotData(:,2),'o-');
    
    plot(distPlotData(idxOfBestPoint(2),1), distPlotData(idxOfBestPoint(2),2), 'og','LineWidth',3);
    legend([ {'all APCA'} {sprintf('Second Best fit %i',nBest2)}]);
    set(gca,'fontsize',20);
    xlabel('# segments (n)');
    ylabel('error');
    title(newname);
    fname=sprintf('%s_error_second.eps',name);
    print(fig5,fname,'-depsc');
    hold off;
    close(fig5);
    %[M,I]=min(distPlotData(:,2));
    
    nBest2
     
    fig6=figure('Visible','off');
    hold all;
    plot(distPlotData(:,1),distPlotData(:,2),'o-');
    
    plot(distPlotData(idxOfBestPoint(3),1), distPlotData(idxOfBestPoint(3),2), 'oc','LineWidth',3);
    legend([ {'all APCA'} {sprintf('Third Best fit %i',nBest3)}]);
    set(gca,'fontsize',20);
    xlabel('# segments (n)');
    ylabel('error');
    title(newname);
    fname=sprintf('%s_error_third.eps',name);
    print(fig6,fname,'-depsc');
    hold off;
    close(fig6);
    %[M,I]=min(distPlotData(:,2));
    
    nBest3
    
    
    fig7=figure('Visible','off');
    hold all;
    plot(distPlotData(:,1),distPlotData(:,2),'o-');
    
    plot(distPlotData(idxOfBestPoint(4),1), distPlotData(idxOfBestPoint(4),2), 'oy','LineWidth',3);
    legend([ {'all APCA'} {sprintf('Fourth Best fit %i',nBest4)}]);
    set(gca,'fontsize',20);
    xlabel('# segments (n)');
    ylabel('error');
    
    title(newname);
    fname=sprintf('%s_error_fourth.eps',name);
    print(fig7,fname,'-depsc');
    hold off;
    close(fig7);
    
    %[M,I]=min(distPlotData(:,2));
    
    nBest4
 %% Plot best segmentation
    fig3=figure('Visible','off');
    hold all;
    plot(distToLine,'o-');
    plot(nBest, distToLine(idxOfBestPoint(1)), 'or','LineWidth',3);
    plot(nBest2, distToLine(idxOfBestPoint(2)), 'og','LineWidth',3);
    plot(nBest3, distToLine(idxOfBestPoint(3)), 'oc','LineWidth',3);
    plot(nBest4, distToLine(idxOfBestPoint(4)), 'oy','LineWidth',3);
    set(gca,'fontsize',16);
    xlabel('# segments (n)');
    ylabel('distance');
    
    title(newname);
    fname=sprintf('%s_%d_%d_%d_%d_segments_distance.eps',name,nBest,nBest2,nBest3,nBest4);
    legend([ {'all APCA segments'} {sprintf('Best fit %i',nBest)} {sprintf('Second best fit %i',nBest2)} {sprintf('Third Best fit %i',nBest3)} {sprintf('Fourth Best fit %i',nBest4)}]);
    print(fig3,fname,'-depsc');
    hold off;
    close(fig3);

    %% Plot dist VS n
    fig4=figure('Visible','off');
    
    hold all; 
    %plot(oMat,'b','LineWidth',1); 
    plot(a*0.1,mMat,'m','LineWidth',1); 
    
    [dAPCA,xAPCA,yAPCA] = aPCAMe(mMat,nBest);
    
    plot(xAPCA*0.1,yAPCA,'r-','LineWidth',3);
    legend([ {'Medfilt data'} {sprintf('APCA Best fit %i',nBest)}],'Location','best');
    set(gca,'fontsize',16);
    xlabel('Genome size (Mb)');
    ylabel('log2(normalized FPKM)');
    mu = mean(pMat(~isinf(pMat)));
    hline = refline([0 mu]);
    %
    
    % build return structure
    cRet.x    = xAPCA;
    cRet.y    = yAPCA;
    cRet.n    = nBest;
    cRet.apca = dAPCA;
    
    cRet.err  = distPlotData;
    title(newname);
    fname=sprintf('%s_%d_segmentation.eps',name,nBest);
    print(fig4,fname,'-depsc');
    hold off;
    close(fig4);
    
    fig8=figure('Visible','off'); 
    
    hold all; 
    %plot(oMat,'b','LineWidth',1); 
    plot(a*0.1,mMat,'m','LineWidth',1); 
    [d2APCA,x2APCA,y2APCA] = aPCAMe(pMat,nBest2);
    plot(x2APCA*0.1,y2APCA,'r-','LineWidth',3);
    legend([ {'Medfilt data'} {sprintf('APCA Second Best fit %i',nBest2)}],'Location','best');
    set(gca,'fontsize',16);
    xlabel('Genome size (Mb)');
    ylabel('log2(normalized FPKM)');
    mu = mean(pMat(~isinf(pMat)));
    hline = refline([0 mu]);
    %
    
    % build return structure
    cRet2.x    = x2APCA;
    cRet2.y    = y2APCA;
    cRet2.n    = nBest2;
    cRet2.apca = d2APCA;
    title(newname);
    cRet2.err  = distPlotData;
    
    fname=sprintf('%s_%d_segmentation.eps',name,nBest2);
    print(fig8,fname,'-depsc');
    % build return structure
    hold off;
    close(fig8);
    
    fig9=figure('Visible','off'); 
    
    hold all; 
    %plot(oMat,'b','LineWidth',1); 
    plot(a*0.1,mMat,'m','LineWidth',1); 
    [d3APCA,x3APCA,y3APCA] = aPCAMe(pMat,nBest3);
    plot(x3APCA*0.1,y3APCA,'r-','LineWidth',3);
    legend([ {'Medfilt data'} {sprintf('APCA Third Best fit %i',nBest3)}],'Location','best');
    %DefaultTextFontSize does not work
    set(gca,'fontsize',16);
    xlabel('Genome size (Mb)');
    ylabel('log2(normalized FPKM)');
    mu = mean(pMat(~isinf(pMat)));
    hline = refline([0 mu]);
    %
    
    % build return structure
    cRet3.x    = x3APCA;
    cRet3.y    = y3APCA;
    cRet3.n    = nBest3;
    cRet3.apca = d3APCA;
    title(newname);
    cRet3.err  = distPlotData;
    
    fname=sprintf('%s_%d_segmentation.eps',name,nBest3);
    print(fig9,fname,'-depsc');
    % build return structure
    hold off;
    close(fig9);
    
    fig10=figure('Visible','off'); hold all; 
    %plot(oMat,'b','LineWidth',1); 
    plot(a*0.1,mMat,'m','LineWidth',1); 
    [d4APCA,x4APCA,y4APCA] = aPCAMe(pMat,nBest4);
    plot(x4APCA*0.1,y4APCA,'r-','LineWidth',3);
    legend([ {'Medfilt data'} {sprintf('APCA Fourth Best fit %i',nBest4)}],'Location','best');
    set(gca,'fontsize',16);
    xlabel('Genome size (Mb)');
    ylabel('log2(normalized FPKM)');
    mu = mean(pMat(~isinf(pMat)));
    hline = refline([0 mu]);
   %
    
    % build return structure
    cRet4.x    = x4APCA;
    cRet4.y    = y4APCA;
    cRet4.n    = nBest4;
    cRet4.apca = d4APCA;
    title(newname);
    cRet4.err  = distPlotData;
    
    fname=sprintf('%s_%d_segmentation.eps',name,nBest4);
    print(fig10,fname,'-depsc');
    % build return structure
    hold off;
    close(fig10);
    
    
    fig11=sfigure(); 
    set(0, 'DefaultFigureVisible', 'off');
    %plot(oMat,'b','LineWidth',1); 
    suptitle(newname);
    set(0, 'currentfigure', fig11);
    subplot(2,2,1);  hold on; 
    plot(a*0.1,mMat,'m','LineWidth',1); 
    xlabel('Genome size (Mb)');
    ylabel('log2(normalized FPKM)');
    plot(xAPCA*0.1,yAPCA,'r-','LineWidth',3);
    
    legend([{'Medfilt data'} {sprintf('APCA Best fit %i',nBest)}],'Location','northoutside');
    set(gca,'FontSize',7);
    mu = mean(pMat(~isinf(pMat)));
    hline = refline([0 mu]);
    
    set(0, 'currentfigure', fig11);
    subplot(2,2,2);  hold on;
    plot(a*0.1,mMat,'m','LineWidth',1); 
    xlabel('Genome size (Mb)');
    ylabel('log2(normalized FPKM)');
    plot(x2APCA*0.1,y2APCA,'r-','LineWidth',3);
    legend([{'Medfilt data'} {sprintf('APCA Second Best fit %i',nBest2)}],'Location','northoutside');
    set(gca,'FontSize',7);
    mu = mean(pMat(~isinf(pMat)));
    hline = refline([0 mu]);
    %
    set(0, 'currentfigure', fig11);
    subplot(2,2,3); hold on; 
    plot(a*0.1,mMat,'m','LineWidth',1); 
    xlabel('Genome size (Mb)');
    ylabel('log2(normalized FPKM)');
    plot(x3APCA*0.1,y3APCA,'r-','LineWidth',3);
    legend([{'Medfilt data'} {sprintf('APCA Third Best fit %i',nBest3)}],'Location','northoutside');
    set(gca,'FontSize',7);
    mu = mean(pMat(~isinf(pMat)));
    hline = refline([0 mu]);
    
    set(0, 'currentfigure', fig11);
    subplot(2,2,4); hold on; 
    plot(a*0.1,mMat,'m','LineWidth',1); 
    xlabel('Genome size (Mb)');
    ylabel('log2(normalized FPKM)');
    plot(x4APCA*0.1,y4APCA,'r-','LineWidth',3);
    legend([{'Medfilt data'} {sprintf('APCA Fourth Best fit %i',nBest4)}],'Location','northoutside');
    set(gca,'FontSize',7);
    
    mu = mean(pMat(~isinf(pMat)));
    hline = refline([0 mu]);
    
    
    fname=sprintf('%s_all_best4_segmentation.eps',name);
    
    print(fig11,fname,'-depsc');

    % build return structure
    hold off;
    close(fig11);
    
    
    %set(fig11,'visible','off');

    %fid=fopen('apca_out.txt','w');
    %fprintf(fid,'%d\t%d\t%d\t%d\n',cRet.apca.lx,cRet.apca.rx,cRet.apca.mc,cRet.apca.y);
end


function [rdAPCA,rxAPCA,ryAPCA] = aPCAMe(apMat,apN)
    
    % Create X,Y matrices
    rxAPCA = zeros(size(apMat));
    ryAPCA = zeros(size(apMat));
    
    % Run APCA
    rdAPCA = apca(apMat',apN);
    
    % Build return X,Y matrices from APCA structure
    for i = 1:length(rdAPCA);
    
        nSt  = rdAPCA(i).lx;
        nSt  = rdAPCA(i).lx;
        nEnd = rdAPCA(i).rx;
    
        rxAPCA(nSt:nEnd) = nSt:nEnd;
        ryAPCA(nSt:nEnd) = ones(1,nEnd-nSt+1).*rdAPCA(i).y;
    end
       
end

function h = sfigure(h)
% SFIGURE  Create figure window (minus annoying focus-theft).
%
% Usage is identical to figure.
%
% Daniel Eaton, 2005
%
% See also figure

if nargin>=1 
	if ishandle(h)
		set(0, 'CurrentFigure', h);
	else
		h = figure(h);
	end
else
	h = figure;
end
end