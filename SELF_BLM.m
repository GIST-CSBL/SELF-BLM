%%%%%%%%%%%%%%%%%%%%% Generate negative labeling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To use k-medodis clustering, the Statistics and Machine Learning Toolbox is used.
% k-medodis clustering 
% N : divid integer  
% target : similiarity matrix of target 
% com : simliarity matrix of compound
% unlabeledY : drug-target interaction matrix
% labeledY : drug- target interaction matrix which has predicted negative interactions
% these original similarity and interaction datasets are avaliable at the
% following URL http://cbio.mines-paristech.fr/~yyamanishi/bipartitelocal/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

target = textread('gpcr_simmat_dg2.txt');
comp = textread('gpcr_simmat_dc2.txt');
unlabeledY = textread('gpcr_admat_dgc2.txt');

N = 2;

target = target(:,2:(size(target,2)));
X = [1:length(target)]';
clus1 = round(size(target,2)/N);
target_idx = kmedoids(X, clus1,'distance',@gene_dist);

comp = comp(:,2:(size(comp,2)));
X = [1:length(comp)]';
clus2 = round(size(comp,2)/N);
com_idx = kmedoids(X,clus2,'distance',@com_dist);

unlabeledY = unlabeledY(:,2:(size(unlabeledY,2)));
unlabeledY = unlabeledY';

labeledY=zeros(size(unlabeledY,1),size(unlabeledY,2));
for i =1:size(unlabeledY,1)
    for j =1 : size(unlabeledY,2)
        if unlabeledY(i,j) ==1
            labeledY(i,j) =1;
        else
            m = 0;
            index = find(target_idx == target_idx(j));
            for k = 1:size(index,1)
                comin = find(unlabeledY(:,index(k)) == 1);
                for q = 1:size(comin,1)
                    if(com_idx(comin(q)) == com_idx(i))
                        labeledY(i,j) = 0;
                        m = 1;
                        break;
                    end
                end
                if(m ==1)
                    break;
                end
            end
            if(m==0)
                labeledY(i,j) = -1;
            end
        end                  
    end
end

labeledY = labeledY';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SELF_BLM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this method is based on bipartite local model (BLM)
% BLM code can find this site : http://cbio.mines-paristech.fr/~yyamanishi/bipartitelocal/
% the libsvm package must be installed

% Purpose : prediction of drug-target interaction. we did Leave-One-Out Cross-Validation

% Output : compPred : Predict drug-target interactions using compound similarity
%          targetPred : Preidction drug-target interactions using target similarity


comp = (comp + comp')/2;

% alpha : thrash hold of selftraining SVM
alpha =1;

% initialize
compLength = length(comp); 
targetLength = length(target);
compPred = zeros(compLength,targetLength);
targetPred = zeros(targetLength,compLength);

%Checking for positive semi-definite
epsilon = .1;
while sum(eig(comp) >= 0) < compLength || isreal(eig(comp))==0 
    comp = comp + epsilon*eye(compLength);
end
while sum(eig(target) >= 0) < targetLength || isreal(eig(target))==0 
    target = target + epsilon*eye(targetLength);
end

for i=1:targetLength
    currentY = labeledY(i,:)';
    for j = 1:compLength

        if sum(currentY(setdiff(1:compLength,j:j)) == 1) > 0 

            trainK = [[1:(compLength-1)]' comp(setdiff(1:compLength,j:j),setdiff(1:compLength,j:j))];
            testK = [(j)' comp(j,setdiff(1:compLength,j:j))];
            kcurrentY = currentY(setdiff(1:compLength,j:j));
            
            if sum(kcurrentY == -1) > 0
               
               trainK1 = trainK(:,2:(size(trainK,2)));
               testK1 = testK(:,2:(size(testK,2))); 
               
               while sum(kcurrentY ==0) >0
                    
                    labeled_index = find(kcurrentY ~= 0);
                    
                    sTrainK = trainK1(labeled_index,labeled_index);
                    sTrainK = [[1:size(labeled_index)]' sTrainK];
                    model = svmtrain(double(kcurrentY(labeled_index)==1),sTrainK,strcat(['-t 4 -c 1 -w1 ','1',' -w-1 ','1']));
                    
                    unlabeled_index = find(kcurrentY == 0);
                    sTestK = trainK1(unlabeled_index,labeled_index);
                    sTestK = [[1:size(unlabeled_index)]' sTestK];
                    [pre_label,acc,dec] = svmpredict(zeros(size(unlabeled_index)),sTestK,model);
                    th = find(abs(dec) > alpha);
                    pre_label1 = pre_label(th);
                    kcurrentY(unlabeled_index(th)) = sign(pre_label1-1/2);
                    if size(unlabeled_index,1) == size(find(kcurrentY == 0),1) | size(find(kcurrentY ==0),1) == 0
                        break;
                    end
               end
                
                labeled_index = find(kcurrentY ~= 0);
                sTrainK = trainK1(labeled_index,labeled_index);
                sTrainK = [[1:size(labeled_index)]' sTrainK];
                model = svmtrain(double(kcurrentY(labeled_index)==1),sTrainK,strcat(['-t 4 -c 1 -w1 ','1',' -w-1 ','1']));
                
                testK1 = testK1(1,labeled_index);
                testK1 = [(j)' testK1];
                [predict_label,accuracy,dec_values] = svmpredict(0,testK1, model);
                firstLabel = double(kcurrentY(labeled_index(1))==1);
                myP1 = dec_values*sign(firstLabel - 1/2);
                compPred(j,i) = myP1;
            else
                
                model = svmtrain(double(kcurrentY==1),trainK,strcat(['-t 4 -c 1 -w1 ','1',' -w-1 ','1']));
                [predict_label,accuracy,dec_values] = svmpredict(0,testK, model);
                firstLabel = double(kcurrentY(1)==1);
                compPred(j,i) = dec_values*sign(firstLabel - 1/2);   
            end
                
        else
            compPred(j,i) = -5; 
        end
    end
end


 for i= 1:compLength

    currentY = labeledY(:,i);

    for j = 1:targetLength
        if sum(currentY(setdiff(1:targetLength,j:j)) == 1) > 0

          trainK = [[1:(targetLength-1)]' target(setdiff(1:targetLength,[j:j]),setdiff(1:targetLength,[j:j]))];
          testK = [(j)' target(j,setdiff(1:targetLength,[j:j]))];
          kcurrentY = currentY(setdiff(1:targetLength,j:j));

          if sum(kcurrentY == -1) > 0

               trainK1 = trainK(:,2:(size(trainK,2)));
               testK1 = testK(:,2:(size(testK,2))); 
              
               while sum(kcurrentY ==0) > 0
                    
                   labeled_index = find(kcurrentY ~= 0);
                    sTrainK = trainK1(labeled_index,labeled_index); 
                    sTrainK = [[1:size(labeled_index)]' sTrainK];
                    model = svmtrain(double(kcurrentY(labeled_index)==1),sTrainK,strcat(['-t 4 -c 1 -w1 ','1',' -w-1 ','1']));
                    
                    unlabeled_index = find(kcurrentY == 0);
                    sTestK = trainK1(unlabeled_index,labeled_index);
                    sTestK = [[1:size(unlabeled_index)]' sTestK];
                    [pre_label,acc,dec] = svmpredict(zeros(size(unlabeled_index)),sTestK,model);
                    th = find(abs(dec) > alpha);
                    pre_label1 = pre_label(th);
                    kcurrentY(unlabeled_index(th)) = sign(pre_label1-1/2);
                    if size(unlabeled_index,1) == size(find(kcurrentY == 0),1) | size(find(kcurrentY ==0),1) == 0
                        break;
                    end
                end

                labeled_index = find(kcurrentY ~= 0);
                trainK1 = trainK(:,2:(size(trainK,2)));
                sTrainK = trainK1(labeled_index,labeled_index);
                sTrainK = [[1:size(labeled_index)]' sTrainK];
                model = svmtrain(double(kcurrentY(labeled_index)==1),sTrainK,strcat(['-t 4 -c 1 -w1 ','1',' -w-1 ','1']));
                
                testK1 = testK(:,2:(size(testK,2)));
                testK1 = testK1(1,labeled_index);
                testK1 = [(j)' testK1];
                [predict_label,accuracy,dec_values] = svmpredict(0,testK1, model);
                firstLabel = double(kcurrentY(labeled_index(1))==1);
                myP1 = dec_values*sign(firstLabel - 1/2);
                targetPred(j,i) = myP1;              
          else
                
                model = svmtrain(double(kcurrentY==1),trainK,strcat(['-t 4 -c 1 -w1 ','1',' -w-1 ','1']));
                [predict_label,accuracy,dec_values] = svmpredict(0,testK, model);
                firstLabel = double(kcurrentY(1)==1);
                targetPred(j,i) = dec_values*sign(firstLabel - 1/2);
          end

        else
            targetPred(j,i) = -5;
        end
    end
 end


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Validation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % caculate the AUC and AUPR of the model
 % v : drug-target interaction matrix
 % we take the largest value between targetPred and compPred
 
v = textread('validation_set_gpcr.txt');
v = v(:,2:size(v,2));

 maxPred = max(targetPred',compPred);

pred =[];
real = [];
for k = 1:compLength
    pred = [pred maxPred(k,:)];
    real = [real v(:,k)'];
end

mRange = (-5):0.001:(5);
rangeLength = length(mRange);

trueP = zeros(1,rangeLength);
falseP = zeros(1,rangeLength);

mIndex = 0;
for moveit = mRange
        mIndex = mIndex + 1;
        trueP(mIndex) = sum(sign(pred + moveit)==1 & real==1);
        falseP(mIndex) = sum(sign(pred + moveit)==1 & real==0);
end


%%% Calculating AUC %%%

xAUC = falseP/max(falseP);  %%%% 1 - Specificity
yAUC = trueP/max(trueP);    %%%% Sensitivity
AUC = trapz(xAUC,yAUC)

%%% Calculating AUPR %%%

xAUPR = yAUC;  %%%% Recall (Sensitivity)
yAUPR = 1- falseP./(falseP + trueP + .000000001*(falseP==0 & trueP==0)); %%%% Precision %%%%
AUPR = trapz(xAUPR,yAUPR)

