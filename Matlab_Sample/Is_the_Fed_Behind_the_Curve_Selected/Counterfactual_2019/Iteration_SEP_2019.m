for i=1:size(Ey,1)
     
    ['Date of SEP forecast: ', num2str((SEP_hor(1,i))), ', Iteration starts in: ', num2str((SEP_hor(1,starting_year)))] % 使用SEP_hor矩阵的第一行中的日期数据，并将其转换为字符串形式，然后输出一个包含日期信息的字符串
    
    if i >= starting_year;

        %Rename median SEP forecasts 和long run forecast的偏差
        EU=Ey_gap(i,1:H,1)';
        EPi=Ey_gap(i,1:H,2)';
        EY=[Ey_gap(i,1:H,1)'; Ey_gap(i,1:H,2)'];
        EYSig=[EySig(i,1:H,1)'; EySig(i,1:H,2)'];
    
        EFFR=Ey(i,1:H,3)';  
        EFFR0=EFFR;
     
        %Create simple correlation structure across horizons
        for ii=1:H
            for jj=1:H
                RHO(ii,jj)=rho^abs(jj-ii);
            end
        end
    
        SIG=diag(EYSig(1:H));
        COV=1*SIG*RHO*SIG;
        EYsim(1:H,:)=mvnrnd(EY(1:H),(COV),TsimEY)';
        SIG=diag(EYSig(H+1:2*H));
        COV=1*SIG*RHO*SIG;
        EYsim(H+1:2*H,:)=mvnrnd(EY(H+1:2*H),(COV),TsimEY)';

        EYlb=prctile(EYsim',2.5)';
        EYub=prctile(EYsim',100-2.5)';
        EPilb=EYlb(H+1:end);  EPiub=EYub(H+1:end);
        EUlb=EYlb(1:H);  EUub=EYub(1:H);
    
        %Weigthing-matrix:
        w=[lambda 0; 0 1];
        W=kron(w,eye(H,H));
        %W=kron(w,diag(0.7.^(0:H-1))); % Consider the discount rate
    
        index_ij=1;
        for ij_EY=1:TsimEY
            for ij_X=1:TsimX
            
                %Two instr
                b=inv([X1sim(:,ij_X) X2sim(:,ij_X)]'*W*[X1sim(:,ij_X) X2sim(:,ij_X)])*[X1sim(:,ij_X) X2sim(:,ij_X)]'*W*EYsim(:,ij_EY);
                B1(ij_EY,ij_X)=-b(1);
                B2(ij_EY,ij_X)=-b(2);
            
            
                %First instr alone
                b1(ij_EY,ij_X)=-inv([X1sim(:,ij_X)]'*W*[X1sim(:,ij_X)])*X1sim(:,ij_X)'*W*EYsim(:,ij_EY);
            
                %Second instr alone
                b2(ij_EY,ij_X)=-inv([X2sim(:,ij_X)]'*W*[X2sim(:,ij_X)])*X2sim(:,ij_X)'*W*EYsim(:,ij_EY);
                

                %Impose ZLB constraint     
                if dum_ZLBcstr==1
                    if (SEP_hor(1,i)>=2020.2 && SEP_hor(1,i)<=2022.2) | (SEP_hor(1,i)>=2008.7 && SEP_hor(1,i)<=2015.8)
                        if B1(ij_EY,ij_X)<=0
                            B1(ij_EY,ij_X)=0;
                            B2(ij_EY,ij_X)= b2(ij_EY,ij_X);
                        end
                    end
                end       
            
                EFFR1sim(:,index_ij)=EFFR0+B1(ij_EY,ij_X)*Xsim(:,3,ij_X)+B2(ij_EY,ij_X)*Xsim(:,6,ij_X);
            
               index_ij=index_ij+1;
            end
        end
    
        temp1= B1(:,:);
        OPP1(:,i)=temp1(:);
    
        temp2= B2(:,:);
        OPP2(:,i)=temp2(:);
    
    
        temp3= b1(:,:);
        OPP_inst1(:,i)=temp3(:); %ffr instrument alone
    
        temp4= b2(:,:);
        OPP_inst2(:,i)=temp4(:); %slope instrument alone
    
    
        %OPP point estimates (no uncertainty)
        %Two instr
        b=inv([X1 X2]'*W*[X1 X2])*[X1 X2]'*W*EY;
        delta1(i)=-b(1);
        delta2(i)=-b(2);
    
        %First instr alone
        b=inv([X1]'*W*[X1])*X1'*W*EY;
        delta_inst1(i)=-b(1);
    
        %Second instr alone
        b=inv([X2]'*W*[X2])*X2'*W*EY;
        delta_inst2(i)=-b(1);

%% Begin Iteration
        B1_flattened = reshape(B1, [], 1);
    
        % Keep all the results of OPP to compute confidence intervals
        for mm = 1:(TsimEY*TsimX)

            J1 = B1_flattened(mm, :);

            deltaY = squeeze(prctile(X1sim,50,2))*J1; % delta EY


            % Fine SEP_hor(1, i) in bTable.yearn's line k
            idx = find(bTable.yearn == SEP_hor(1, i));

            if ~isempty(idx)
                % If the number of the lines left is no more than 20, cope with the lines left
                rows_to_process = min(height(bTable) - idx + 1, H);

                % Add the first (rows_to_process) lines of delta_Y to bTable.UR's line k to line (k+rows_to_process-1)
               
                % The upper code is point estimation
                % bTable.UR(idx:idx+rows_to_process-1) = bTable.UR(idx:idx+rows_to_process-1) + deltaY(1:rows_to_process);
                bTable_output(idx:idx+rows_to_process-1, 4, mm) = bTable_output(idx:idx+rows_to_process-1, 4, mm) + deltaY(1:rows_to_process);

                % Add the last (rows_to_process) lines of delta_Y to bTable.PI's line k to line (k+rows_to_process-1)

                % Th eupper line is point estimation
                % bTable.PI(idx:idx+rows_to_process-1) = bTable.PI(idx:idx+rows_to_process-1) + deltaY(21:20+rows_to_process);
                bTable_output(idx:idx+rows_to_process-1, 5, mm) = bTable_output(idx:idx+rows_to_process-1, 5, mm) + deltaY(H+1:H+rows_to_process);
                
                % FFR's value after OPP

                %bTable.FFR(idx:idx+rows_to_process-1) = bTable.FFR(idx:idx+rows_to_process-1) + J1;
                bTable_output(idx:idx+rows_to_process-1, 6, mm) = bTable_output(idx:idx+rows_to_process-1, 6, mm) + X3(1:rows_to_process)*J1;

                % impose ZLB constraint
                rows_to_check = find((bTable_output(:, 1, mm) >= 2020.25 & bTable_output(:, 1, mm) <= 2022.25) | (bTable_output(:, 1, mm) >= 2008 & bTable_output(:, 1, mm) <= 2015.75));
                for row = rows_to_check'
                    if bTable_output(row, 6, mm) <= 0
                        bTable_output(row, 6, mm) = 0;
                    end
                end

            else
                
                % Error sign
                disp('No corresponding row found in bTable.yearn for SEP_hor(1, i).');
            end
        end

        % Compress bTable_ouput to bTable for convenience of revising
        % backcast data
        bTable = median(bTable_output, 3);
        colNames = {'yearn','year','quarter','UR','PI','FFR'};
        bTable = array2table(bTable,'VariableNames',colNames);


        %-------------Revise forecast data (t0, t1, t2, t3)---------------%
        if (SEP_hor(1, i) >= 2020.25 & SEP_hor(1, i) <= 2022.25) | (SEP_hor(1, i) >= 2008 & SEP_hor(1, i) <= 2015.75)
            if J1 <= 0
                J1 = 0;
            end
        end

        IRfore = squeeze(prctile(X1sim,50,2));
        forecast_adj = IRfore * J1;

        FFRforecast_adj = squeeze(prctile(squeeze(IR(:,3,INDrnd)),50,2)) * J1;


        % Coefficient of forecast adjustment
        URforecast_adj = 0.18*forecast_adj(1:H-1)';
        PIforecast_adj = 0.13*forecast_adj(H+1:2*H-1)';
        FFRforecast_adj = 0.74*FFRforecast_adj(1:H-1)';


        if i <= numel(SEP_hor(1, :)) - 1
            for i_fore = i: min( numel(SEP_hor(1, :)), i+H-1)
                num_to_drop = round((SEP_hor(1, i_fore) - SEP_hor(1, i)) * 4);
                if num_to_drop <= H-1
                    URforecast_adj_temp = URforecast_adj(num_to_drop + 1 : end);
                    URforecast_adj_temp = horzcat(URforecast_adj_temp, zeros(1, num_to_drop + 1));
                    Ey_gap(i_fore, :, 1) = Ey_gap(i_fore, :, 1) + URforecast_adj_temp;
                else
                end
            end

            for i_fore = i : min( numel(SEP_hor(1, :)), i+H-1)
                num_to_drop = round((SEP_hor(1, i_fore) - SEP_hor(1, i)) * 4);
                if num_to_drop <= H-1
                    PIforecast_adj_temp = PIforecast_adj(num_to_drop + 1 : end);
                    PIforecast_adj_temp = horzcat(PIforecast_adj_temp, zeros(1, num_to_drop + 1));
                    Ey_gap(i_fore, :, 2) = Ey_gap(i_fore, :, 2) + PIforecast_adj_temp;
                else
                end
            end

            for i_fore = i : min( numel(SEP_hor(1, :)), i+H-1)
                if numel(SEP_hor(1, :)) - i_fore + 1 <= 33
                    num_to_drop = round((SEP_hor(1, i_fore) - SEP_hor(1, i)) * 4);
                    if num_to_drop <= H-1
                        FFRforecast_adj_temp = FFRforecast_adj(num_to_drop + 1 : end);
                        FFRforecast_adj_temp = horzcat(FFRforecast_adj_temp, zeros(1, num_to_drop + 1));
                        Ey_gap(i_fore, :, 3) = Ey_gap(i_fore, :, 3) + FFRforecast_adj_temp;
                    else
                    end
                else
                end
            end
        end
        %-----------------------------------------------------------------%

        %-----Fill in the backcast data to simulate forecast path---------%
        for i = length(name)
            % UR
            if i == 1
                RTlr=sTable_UR.median_lr;
                RTlr(isnan(sTable_UR.median_lr))=sTable_UR.RealTime_LR(isnan(sTable_UR.median_lr));
                dataEu=[ sTable_UR.median_t sTable_UR.median_t1 sTable_UR.median_t2 sTable_UR.median_t3 RTlr];

                sTable_UR.backcast = NaN(height(sTable_UR), 1);

                for j = 1:height(sTable_UR)
                    current_year = sTable_UR.year(j);
                    idx = find(bTable.year == current_year - 1 & bTable.quarter == 4);
                    if ~isempty(idx)
                        sTable_UR.backcast(j) = bTable.UR(idx);
                    end
                end

                dataEub=sTable_UR.backcast; 
                dataEubActual=sTable_UR.backcast;

                Eu=NaN(length(dataEu),H);
                for kk=1:length(dataEu)
        
                    if month(kk)==1, j=3; jb=0; end;
                    if month(kk)==3 | month(kk)==4, j=3; jb=0; end;
                    if month(kk)==6, j=2; jb=1; end;
                    if month(kk)==9, j=1; jb=2; end;
                    if month(kk)==12 | month(kk)==11, j=0; jb=3; end;
        
                    for ii=1:4
                        Eu(kk,j+4*(ii-1)+1)=dataEu(kk,ii);
                    end
        
                    %fill the long-run value for year 5
                    Eu(kk,end)=dataEu(kk,end);
        
                    if isnan(dataEub(kk))~=1
                        ExtEu{kk}=[dataEub(kk) NaN(1,jb) Eu(kk,:)];
                    else
                        ExtEu{kk}=[dataEubActual(kk) NaN(1,jb) Eu(kk,:)];
                    end
        
                    temp1=naninterp(ExtEu{kk});
        
                    %remove backcast to start at t=0
                    Eui(kk,:)=temp1(2+jb:end);
        
        
                    %Forecast for the mandates in deviations from their targets (the real-time estimates of the long-run values)
                    Eui_gap(kk,:)=Eui(kk,:)-dataEu(kk,end);
        
                    %date and horizon associated with forecast
                    SEP_hor(:,kk)=(year(kk)+((quarter(kk)-1)/4):1/4:year(kk)+((quarter(kk)-1)/4)+(H-1)/4)';
        
                end
    
                %store output
                Ey(:,:,i)=Eui;
                Ey_gap(:,:,i)=Eui_gap;
            else
                % PI
                RTlr=sTable_PI.median_lr;
                RTlr(isnan(sTable_PI.median_lr))=sTable_PI.RealTime_LR(isnan(sTable_PI.median_lr));
                dataEu=[ sTable_PI.median_t sTable_PI.median_t1 sTable_PI.median_t2 sTable_PI.median_t3 RTlr];
                
                sTable_PI.backcast = NaN(height(sTable_PI), 1);

                for j = 1:height(sTable_PI)
                    current_year = sTable_PI.year(j);
                    idx = find(bTable.year == current_year - 1 & bTable.quarter == 4);
                    if ~isempty(idx)
                        sTable_PI.backcast(j) = bTable.PI(idx);
                    end
                end

                dataEub=sTable_PI.backcast; 
                dataEubActual=sTable_PI.backcast;

                Eu=NaN(length(dataEu),H);
 
                for kk=1:length(dataEu)
        
                    if month(kk)==1, j=3; jb=0; end;
                    if month(kk)==3 | month(kk)==4, j=3; jb=0; end;
                    if month(kk)==6, j=2; jb=1; end;
                    if month(kk)==9, j=1; jb=2; end;
                    if month(kk)==12 | month(kk)==11, j=0; jb=3; end;
        
                    for ii=1:4
                        Eu(kk,j+4*(ii-1)+1)=dataEu(kk,ii);
                    end
        
                    %fill the long-run value for year 5
                    Eu(kk,end)=dataEu(kk,end);
        
                    if isnan(dataEub(kk))~=1
                        ExtEu{kk}=[dataEub(kk) NaN(1,jb) Eu(kk,:)];
                    else
                        ExtEu{kk}=[dataEubActual(kk) NaN(1,jb) Eu(kk,:)];
                    end
        
                    temp1=naninterp(ExtEu{kk});
        
                    %remove backcast to start at t=0
                    Eui(kk,:)=temp1(2+jb:end);
        
        
                    %Forecast for the mandates in deviations from their targets (the real-time estimates of the long-run values)
                    Eui_gap(kk,:)=Eui(kk,:)-dataEu(kk,end);
        
                    %date and horizon associated with forecast
                    SEP_hor(:,kk)=(year(kk)+((quarter(kk)-1)/4):1/4:year(kk)+((quarter(kk)-1)/4)+(H-1)/4)';
        
                end
    
                %store output
                Ey(:,:,i)=Eui;
                Ey_gap(:,:,i)=Eui_gap;

            end
        end
        
        % FFR
        %As long-run value, use SEP and if not available, use real-time long-run estimate
        RTlr=sTable_tp.median_lr;
        RTlr(isnan(sTable_tp.median_lr))=sTable_tp.RealTime_LR(isnan(sTable_tp.median_lr));

        %forecast data are as follows: last Q4 (backcast data), this year's Q4, next year's Q4, t+2yr, t+3yr and long run (t+5yr say)
        dataEu=[ sTable_tp.median_t sTable_tp.median_t1 sTable_tp.median_t2 sTable_tp.median_t3 RTlr];

        sTable_tp.backcast = NaN(height(sTable_tp), 1);

        for j = 1:height(sTable_tp)
            current_year = sTable_tp.year(j);
            idx = find(bTable.year == current_year - 1 & bTable.quarter == 4);
            if ~isempty(idx)
                sTable_tp.backcast(j) = bTable.FFR(idx);
            end
        end

        idx = sTable_tp.year < 2015 | (sTable_tp.year == 2015 & sTable_tp.month <= 6);
        sTable_tp.backcast(idx) = NaN;

        dataEub=[]; 
        dataEubActual=sTable_tp.backcast;

        dataEU_now=sTable_tp.median_t0;

        Eu=NaN(length(dataEu),H);
        Eu_lb=NaN(length(dataEu),H);Eu_ub=NaN(length(dataEu),H);
        Eu_lb2=NaN(length(dataEu),H);Eu_ub2=NaN(length(dataEu),H);
        EuSig=NaN(length(dataEu),H);

        year_tp=sTable_tp.year;
        yearn_tp=sTable_tp.yearn;
        month_tp=sTable_tp.month;

        Eui=[];

        for kk=1:length(dataEu)
    
            %start filling the forecast at different horizon depending on month of FOMC (using "j"), and similarly to fill the backcast ("jb")
            if month_tp(kk)==1, j=3; jb=0; end;
            if month_tp(kk)==3 | month_tp(kk)==4, j=3; jb=0; end;
            if month_tp(kk)==6, j=2; jb=1; end;
            if month_tp(kk)==9, j=1; jb=2; end;
            if month_tp(kk)==12 | month_tp(kk)==11, j=0; jb=3; end;
    
            for ii=1:4
                Eu(kk,j+4*(ii-1)+1)=dataEu(kk,ii);
            end
    
            %fill the long-run value for year 5
            Eu(kk,end)=dataEu(kk,end);
    
            Eui(kk,:)=Eu(kk,:);
    
            %now add current reading (since it could be different from linear interpolation)
            Eui(kk,1)=dataEU_now(kk);
    
            Eui(kk,:)=naninterp( Eui(kk,:));
    
        end

        %store output
        Ey(:,:,3)=Eui;
        Ey_gap(:,:,3)=Eui;
    


        %-----------------End of forecast path adjustment-----------------%
    else
        % Compute the OPP before the starting period, without iteration
        % Rename median SEP forecasts
        EU=Ey_gap(i,1:H,1)';
        EPi=Ey_gap(i,1:H,2)';
        EY=[Ey_gap(i,1:H,1)'; Ey_gap(i,1:H,2)'];
        EYSig=[EySig(i,1:H,1)'; EySig(i,1:H,2)'];
    
        EFFR=Ey(i,1:H,3)';  
        EFFR0=EFFR;
     
        % Create simple correlation structure across horizons
        for ii=1:H
            for jj=1:H
                RHO(ii,jj)=rho^abs(jj-ii);
            end
        end
    
        SIG=diag(EYSig(1:H));
        COV=1*SIG*RHO*SIG;
        EYsim(1:H,:)=mvnrnd(EY(1:H),(COV),TsimEY)';
        SIG=diag(EYSig(H+1:2*H));
        COV=1*SIG*RHO*SIG;
        EYsim(H+1:2*H,:)=mvnrnd(EY(H+1:2*H),(COV),TsimEY)';

        EYlb=prctile(EYsim',2.5)';
        EYub=prctile(EYsim',100-2.5)';
        EPilb=EYlb(H+1:end);  EPiub=EYub(H+1:end);
        EUlb=EYlb(1:H);  EUub=EYub(1:H);
    
        % Weigthing-matrix:
        w=[lambda 0; 0 1];
        W=kron(w,eye(H,H));
        % W=kron(w,diag(0.7.^(0:H-1))); % Consider discount rate
    
        index_ij=1;
        for ij_EY=1:TsimEY
            for ij_X=1:TsimX
            
                %Two instr
                b=inv([X1sim(:,ij_X) X2sim(:,ij_X)]'*W*[X1sim(:,ij_X) X2sim(:,ij_X)])*[X1sim(:,ij_X) X2sim(:,ij_X)]'*W*EYsim(:,ij_EY);
                B1(ij_EY,ij_X)=-b(1);
                B2(ij_EY,ij_X)=-b(2);
            
            
                %First instr alone
                b1(ij_EY,ij_X)=-inv([X1sim(:,ij_X)]'*W*[X1sim(:,ij_X)])*X1sim(:,ij_X)'*W*EYsim(:,ij_EY);
            
                %Second instr alone
                b2(ij_EY,ij_X)=-inv([X2sim(:,ij_X)]'*W*[X2sim(:,ij_X)])*X2sim(:,ij_X)'*W*EYsim(:,ij_EY);
                
                %Impose ZLB constraint     
                if dum_ZLBcstr==1
                    if (SEP_hor(1,i)>=2020.2 && SEP_hor(1,i)<=2022.2) | (SEP_hor(1,i)>=2008.7 && SEP_hor(1,i)<=2015.8)
                        if B1(ij_EY,ij_X)<=0
                            B1(ij_EY,ij_X)=0;
                            B2(ij_EY,ij_X)= b2(ij_EY,ij_X);
                        end
                    end
                end       
            
                EFFR1sim(:,index_ij)=EFFR0+B1(ij_EY,ij_X)*Xsim(:,3,ij_X)+B2(ij_EY,ij_X)*Xsim(:,6,ij_X);
            
               index_ij=index_ij+1;
            end
        end
    
        temp1= B1(:,:);
        OPP1(:,i)=temp1(:);
    
        temp2= B2(:,:);
        OPP2(:,i)=temp2(:);
    
    
        temp3= b1(:,:);
        OPP_inst1(:,i)=temp3(:); %ffr instrument alone
    
        temp4= b2(:,:);
        OPP_inst2(:,i)=temp4(:); %slope instrument alone
    
    
        %OPP point estimates (no uncertainty)
        %Two instr
        b=inv([X1 X2]'*W*[X1 X2])*[X1 X2]'*W*EY;
        delta1(i)=-b(1);
        delta2(i)=-b(2);
    
        %First instr alone
        b=inv([X1]'*W*[X1])*X1'*W*EY;
        delta_inst1(i)=-b(1);
    
        %Second instr alone
        b=inv([X2]'*W*[X2])*X2'*W*EY;
        delta_inst2(i)=-b(1);
    end
end

if saveres==1
    save FedHistory_lambda1_BVAR_ZLB
end

%% Report results:
bTable_output_SEP = bTable_output;
if saver == 1
    save('bTable_output_SEP_2019.mat', 'bTable_output_SEP', 'yearn','OPP1','yearn_tp');
    save('Welfare.mat','Ey_gap','Ey_gap_raw', 'yearn');
end
Counterfactual_inflation_rate_plot_SEP_2019