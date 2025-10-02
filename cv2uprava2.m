% PGA vs SGA porovnanie 
clc; clear; close all;

numgen=1500;      % number of generations
lpop=50;          % number of chromosomes in population
lstring=10;       % number of genes in a chromosome
M=500;            % maximum of the search space

runs=5;      

bestPGA = zeros(1, runs);
bestSGA = zeros(1, runs);
allPGA = zeros(runs, numgen);   
allSGA = zeros(runs, numgen);
avgPGA = zeros(1, numgen);
avgSGA = zeros(1, numgen);

figure(1); clf; hold on; grid on;
title('PGA vs SGA - priemer a najlepsie behy');
xlabel('Generacia'); ylabel('Fitness');

for r = 1:runs

    Space = [ones(1,lstring)*(-M); ones(1,lstring)*M];
    Delta = Space(2,:)/100;

    % PGA 
    gcpobj = parpool(5,"AutoAddClientPath",true,"SpmdEnabled",true);
    addAttachedFiles(gcpobj,"genetic"); 

    perIsland = zeros(5, numgen);

    migInterval = 50;    
    migSize = 5;          

    spmd
        Pop = genrpop(lpop, Space);
        for gen = 1:numgen

            Fit = eggholder(Pop);

            perIsland(spmdIndex, gen) = min(Fit);

            Best=selbest(Pop,Fit,[3,0]);
            Old=selrand(Pop,Fit,17);

            Work1=selsus(Pop,Fit,10);
            Work2=selsus(Pop,Fit,20);
            Work1=crossov(Work1,1,0);

            Work2=mutx(Work2,0.2,Space);
            Work2=muta(Work2,0.2,Delta,Space);

            Pop=[Best;Old;Work1;Work2];

            if mod(gen, migInterval) == 0
                
                migrants = selbest(Pop, Fit, [migSize,0]);
                
                nextIsland = mod(spmdIndex, numlabs) + 1;  
                spmdSend(migrants, nextIsland);      
               
                incoming = spmdReceive(mod(spmdIndex-2,numlabs)+1);
                
                [~, worstIdx] = sort(Fit, 'descend');
                Pop(worstIdx(1:migSize), :) = incoming;
            end
        end
    end
    pData = cat(1, perIsland{:});  
    mergedPGA = min(pData, [], 1);
    allPGA(r,:) = mergedPGA;
    bestPGA(r) = min(mergedPGA);
    avgPGA = avgPGA + mergedPGA;
    delete(gcpobj);

    % SGA 
    sgaPopSize = 5*lpop;
    PopSGA = genrpop(sgaPopSize, Space);
    sgaEvol = zeros(1,numgen);

    for gen = 1:numgen

        FitSGA = eggholder(PopSGA);

        sgaEvol(gen) = min(FitSGA);

        Best = selbest(PopSGA, FitSGA, [5,0]);
        Old  = selrand(PopSGA, FitSGA, sgaPopSize - 30);

        Work1 = selsus(PopSGA, FitSGA, 5);
        Work2 = selsus(PopSGA, FitSGA, 20);
        Work1 = crossov(Work1,1,0);

        Work2 = mutx(Work2,0.15,Space);
        Work2 = muta(Work2,0.15,Delta,Space);

        PopSGA = [Best; Old; Work1; Work2];
        
        if size(PopSGA,1) > sgaPopSize
            PopSGA = PopSGA(1:sgaPopSize,:);
        elseif size(PopSGA,1) < sgaPopSize
            PopSGA = [PopSGA; genrpop(sgaPopSize - size(PopSGA,1), Space)];
        end
    end
    allSGA(r,:) = sgaEvol;
    bestSGA(r) = min(sgaEvol);
    avgSGA = avgSGA + sgaEvol;

    plot(mergedPGA, 'k', 'LineWidth', 0.5);
    plot(sgaEvol, 'k', 'LineWidth', 0.5);
end

% priemery
avgPGA = avgPGA / runs;
avgSGA = avgSGA / runs;

% najlepsie behy
[~, idxBestP] = min(bestPGA);
[~, idxBestS] = min(bestSGA);

% vykreslenie 
hPGAavg = plot(avgPGA, 'Color',[0 1 1], 'LineWidth',4);        
hPGAmin = plot(allPGA(idxBestP,:), 'Color',[0.5 0 0.5], 'LineWidth',2); 
hSGAavg = plot(avgSGA, 'r', 'LineWidth',2);                    
hSGAmin = plot(allSGA(idxBestS,:), 'g', 'LineWidth',2);        

legend([hPGAavg hPGAmin hSGAavg hSGAmin], ...
       {'Priemer PGA','Najlepsi PGA','Priemer SGA','Najlepsi SGA'}, ...
       'Location','northeastoutside');

fprintf('\nPGA: mean of best = %.6g, global best = %.6g\n', mean(bestPGA), min(bestPGA));
fprintf('SGA: mean of best = %.6g, global best = %.6g\n', mean(bestSGA), min(bestSGA));
