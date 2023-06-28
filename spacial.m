function result = spacial(filename, plist,humanlist,mosquitolist,EIRS,prevalences,plotgif,sample,findhuman,findmosquito)
% filename spacial.m
% Main Function
% Arguments:
%   filename: plot save as ./filename.gif
%   plist: parameters saved as struct
%   parameters in plist:
%       drive_conversion: probability for drive-gene to convert wild-type
%       to drived
%       drive_fitness: fitness for dr-homozygote, which heterozygotes with
%       sqrt(df)
%       release_rate: 释放量/环境容纳量
%       germline_resistance_forming: 产生resistance配子的生成率
%       dd_mothertal_inheritance: 基因驱动效应通过母系遗传给下一代的概率
    result = 0x7FFFFFFF;
    storename=filename;
    storename = sprintf('%s_%.3f_%.3f_%.3f_%.3f_%.3f_%.3f_%.3f', storename, plist.drive_conversion, plist.drive_fitness, plist.release_rate, plist.germline_resistance_forming, plist.dd_mothertal_inheritance, plist.log_immunity_speed, plist.reducehtm);
    %drive_fitness,drive_conversion,speed
    
    %% Default parameter
    %times and spatial
    time=800;
    drive_release_time=10;
    release_rate=plist.release_rate;
    n=100;
    
    if sample
        samplepoint1=1;
        samplepoint2=n;
        t1=1000000;
        t2=1000000;
    end
    speed=1;
    if findhuman
        n=length(EIRS);
        speed=0;
    end
    if findmosquito
        n=length(prevalences);
        speed=0;
    end
    
    % competition
    %capacity=reshape((1:n)*0.1,n,1);
    capacity=mosquitolist.capacity;
    old_larva_competition_factor=5;
    Low_density_growth_rate=9;
    num_eggs=16;
    reproduction_rate=mosquitolist.reproduction_rate;
    %reproduction_rate=reshape(0.8*(1:n)/n,n,1);
    %capacity=1./reproduction_rate./reproduction_rate./reproduction_rate;
    %malaria
    
    immunity_speed=10^plist.log_immunity_speed;
    htom=reproduction_rate*mosquitolist.htom; %rate that healthy become incubation per human rate
    develop_rate=5/8; %rate that incubation become patient
    %develop_rate=reshape((1:n)/n,n,1);
    x=1:n;
    
    population=1;
    humanlist.mtoh=humanlist.mtoh*reproduction_rate/population;

    %humanlist.human_recovery=0.05;%0.1-0.03
    human_resistance=ones(n,1)*2.73;%initial resistance
    humanlist.immunity_gain_rate=humanlist.immunity_gain_rate*immunity_speed;
    humanlist.immunity_losing_rate=humanlist.immunity_losing_rate*immunity_speed;
    %humanlist.b1=0.5;
    %humanlist.shape=2.155;
    
    %malaria drive
    reducedevelop1=0;
    reducedevelop2=0;
    reducehtm1=0.55;
    reducehtm2=0.65;
    reducemth1=0;
    reducemth2=0;
    
    %resistance forms before drive conversion
    drive_conversion=plist.drive_conversion*(1-plist.germline_resistance_forming);
    drive_fitness=plist.drive_fitness;
    germline_resistance_forming=plist.germline_resistance_forming;
    partial_HDR_rate=0.0001;
    NHEJ_r1_forming_rate=0.000001;
    dd_mothertal_inheritance=plist.dd_mothertal_inheritance;
    
    %drive type
    recessivelethal=true;
    haplolethal=false;
    suppression=false;
    

    
    %% mosquito initialize
    
    expected_old=capacity/3;
    expected_larva=capacity*num_eggs*reproduction_rate;
    
    
    %dimension1:spatial
    e=1;
    while normcdf(-e+0.5,0,speed)>0.05/e
        e=e+1;
    end
    
    fussion=1:e;
    for i =1:e
        fussion(i)=normcdf(0.5-i,0,speed);
    end
    
    for i =1:e-1
        fussion(i)=(fussion(i)-fussion(i+1));
    end
    
    %dimension2:genotype
    % wtwt wtdr wtr1 wtr2 drdr drr1 drr2 r1r1 r1r2 r2r2
    
    %dimension3:malaria
    
    %dimension4:time
    mosquitos=zeros(n,10,3,8);
    alive=[6/7,5/6,4/5,3/4,2/3,1/2];
    alivematrix=ones(n,10,3).*reshape(alive,1,1,1,6);
    
    %drive_conversion+resistance_forming should be less than 1;
    
    drive_hete_fitness=sqrt(drive_fitness);
    %for convert allele 
    germline_r1_forming=germline_resistance_forming*NHEJ_r1_forming_rate+partial_HDR_rate;
    germline_r2_forming=germline_resistance_forming*(1-NHEJ_r1_forming_rate);
    male_drive_conversion=drive_conversion;
    female_drive_conversion=drive_conversion;
    
    %for convert another gamete
    
    dn_mothertal_inheritance=dd_mothertal_inheritance*1.87*0.5;
    dn_embryo_normal=1-dn_mothertal_inheritance;
    dn_embryo_r1_forming=dn_mothertal_inheritance*NHEJ_r1_forming_rate;
    dn_embryo_r2_forming=dn_mothertal_inheritance*(1-NHEJ_r1_forming_rate);
    dd_embryo_normal=1-dd_mothertal_inheritance;
    dd_embryo_r1_forming=dd_mothertal_inheritance*NHEJ_r1_forming_rate;
    dd_embryo_r2_forming=dd_mothertal_inheritance*(1-NHEJ_r1_forming_rate);
    
    %% female matrix
    % wtwt wtdr wtr1 wtr2 drdr drr1 drr2 r1r1 r1r2 r2r2
    femalematrix=zeros(1,10,12);
    
    %wtwt
    femalematrix(1,1,1)=1;
    
    %wtdr
    femalematrix(1,2,2)=0.5+0.5*female_drive_conversion;
    femalematrix(1,2,3)=0.5*germline_r1_forming;
    femalematrix(1,2,4)=0.5*germline_r2_forming;
    femalematrix(1,2,1)=1-femalematrix(1,2,2)-femalematrix(1,2,3)-femalematrix(1,2,4);
    femalematrix(1,2,5:8)=femalematrix(1,2,1:4)*dn_embryo_r1_forming;
    femalematrix(1,2,9:12)=femalematrix(1,2,1:4)*dn_embryo_r2_forming;
    femalematrix(1,2,1:4)=femalematrix(1,2,1:4)*dn_embryo_normal;
    
    %wtr1
    femalematrix(1,3,1)=0.5;
    femalematrix(1,3,3)=0.5;
    
    %wtr2
    femalematrix(1,4,1)=0.5;
    femalematrix(1,4,4)=0.5;
    
    %drdr
    if ~suppression
        femalematrix(1,5,2)=dd_embryo_normal;
        femalematrix(1,5,6)=dd_embryo_r1_forming;
        femalematrix(1,5,10)=dd_embryo_r2_forming;
    end
    
    
    
    
    %drr1
    femalematrix(1,6,2)=0.5;
    femalematrix(1,6,3)=0.5;
    femalematrix(1,6,5:8)=femalematrix(1,6,1:4)*dn_embryo_r1_forming;
    femalematrix(1,6,9:12)=femalematrix(1,6,1:4)*dn_embryo_r2_forming;
    femalematrix(1,6,1:4)=femalematrix(1,6,1:4)*dn_embryo_normal;
    
    
    %drr2
    if ~suppression
        femalematrix(1,7,2)=0.5;
        femalematrix(1,7,4)=0.5;
        femalematrix(1,7,5:8)=femalematrix(1,7,1:4)*dn_embryo_r1_forming;
        femalematrix(1,7,9:12)=femalematrix(1,7,1:4)*dn_embryo_r2_forming;
        femalematrix(1,7,1:4)=femalematrix(1,7,1:4)*dn_embryo_normal;
    end
    %r1r1
    femalematrix(1,8,3)=1;
    
    %r1r2
    femalematrix(1,9,3)=0.5;
    femalematrix(1,9,4)=0.5;
    
    %r2r2
    femalematrix(1,10,4)=1;
    
    %% malematrix
    malematrix=zeros(1,10,4);
    
    %wtwt
    malematrix(1,1,1)=1;
    
    %wtdr
    malematrix(1,2,2)=0.5+0.5*male_drive_conversion;
    malematrix(1,2,3)=0.5*germline_r1_forming;
    malematrix(1,2,4)=0.5*germline_r2_forming;
    malematrix(1,2,1)=1-malematrix(1,2,2)-malematrix(1,2,3)-malematrix(1,2,4);
    
    %wtr1
    malematrix(1,3,1)=0.5;
    malematrix(1,3,3)=0.5;
    
    %wtr2
    malematrix(1,4,1)=0.5;
    malematrix(1,4,4)=0.5;
    
    %drdr
    malematrix(1,5,2)=1;
    
    %drr1
    malematrix(1,6,2)=0.5;
    malematrix(1,6,3)=0.5;
    
    %drr2
    malematrix(1,7,2)=0.5;
    malematrix(1,7,4)=0.5;
    
    %r1r1
    malematrix(1,8,3)=1;
    
    %r1r2
    malematrix(1,9,3)=0.5;
    malematrix(1,9,4)=0.5;
    
    %r2r2
    malematrix(1,10,4)=1;
    
    
    
    %% stable loops
    
    
    for i = 1:10
        %aging
        mosquitos(:,:,:,3:8)=mosquitos(:,:,:,2:7).*alivematrix;
        mosquitos(:,:,:,2)=mosquitos(:,:,:,1);
        %genotype initialize   
        new_larva=zeros(n,10,3);
        new_larva(:,1,1)=expected_old;
    
        mosquitos(:,:,:,1)=new_larva;
        %infection
        mosquitos(:,:,3,3:8)=mosquitos(:,:,2,3:8)+mosquitos(:,:,3,3:8);
        mosquitos(:,:,2,3:8)=mosquitos(:,:,1,3:8).*0.4.*htom;
        mosquitos(:,:,3,3:8)=mosquitos(:,:,2,3:8)*develop_rate+mosquitos(:,:,3,3:8);
        mosquitos(:,:,2,3:8)=mosquitos(:,:,2,3:8).*(1-develop_rate);
        mosquitos(:,:,1,3:8)=mosquitos(:,:,1,3:8).*(1-0.4.*htom);
    end
    
    %% human initialize
    
    %dimension1:spatial
    %dimension2:malaria healthy-incubation1-incubation2-malaria
    
    humans=ones(n,4);
    humans(:,1)=humans(:,1)*0.44;
    humans(:,2)=humans(:,2)*0.05;
    humans(:,3)=humans(:,3)*0.05;
    humans(:,4)=humans(:,4)*0.46;
    
    decrease_time=5;
    decrease_speed=decrease_time/2;
    %datacollecting_cycle:
    
    
    mosquito_malarias=zeros(n,time);
    human_malarias=zeros(n,time);
    human_immunitys=zeros(n,time);
    mosquito_genotypes=zeros(n,10,time);
    
    %% main loop
    if plotgif
        f1=figure;
    end
    max=(15).*ones(n,1);
    hmax=ones(n,1)*1.2;
    hmin=zeros(n,1);
        
    for i = 1:time
        if i==drive_release_time
            mosquitos(1:10,2,1,3)=mosquitos(1:10,2,1,3)+release_rate*capacity;
        end
        %aging
        mosquitos(:,:,:,3:8)=mosquitos(:,:,:,2:7).*alivematrix;
        mosquitos(:,:,:,2)=mosquitos(:,:,:,1);
    %%
        % generate new_larvas
        genotypes=sum(sum(mosquitos(:,:,:,3:8),3),4);
        mosquito_genotypes(:,:,i)=genotypes;
        all_numbers=sum(genotypes,2);
        new_larva=zeros(n,10,3);
        egg_cell=reshape(sum(genotypes.*femalematrix,2),n,12);
        sperm=reshape(sum(genotypes./(all_numbers.*ones(1,10)).*malematrix,2),n,1,4);
        eggs=sperm.*egg_cell;
    
        % wtwt wtdr wtr1 wtr2 drdr drr1 drr2 r1r1 r1r2 r2r2
        new_larva(:,3,1)=eggs(:,5,1); %wtr1
        new_larva(:,6,1)=eggs(:,6,1); %drr1
        new_larva(:,8,1)=eggs(:,7,1); %r1r1
        new_larva(:,9,1)=eggs(:,8,1)+eggs(:,11,1); %r1r2
        new_larva(:,4,1)=eggs(:,9,1);
        new_larva(:,7,1)=eggs(:,10,1);
        new_larva(:,10,1)=eggs(:,12,1);
    
        eggs(:,1:4,2:4)=eggs(:,1:4,2:4)+eggs(:,5:8,2:4)+eggs(:,9:12,2:4);
        new_larva(:,1:4,1)=new_larva(:,1:4,1)+reshape(eggs(:,1,1:4),n,4,1);
        new_larva(:,5:7,1)=new_larva(:,5:7,1)+reshape(eggs(:,2,2:4),n,3,1);
        new_larva(:,8:9,1)=new_larva(:,8:9,1)+reshape(eggs(:,3,3:4),n,2,1);
        new_larva(:,10,1)=new_larva(:,10,1)+eggs(:,1,4);
        new_larva(:,2:4,1)=new_larva(:,2:4,1)+eggs(:,2:4,1);
        new_larva(:,6:7,1)=new_larva(:,6:7,1)+eggs(:,3:4,2);
        new_larva(:,9,1)=new_larva(:,9,1)+eggs(:,4,3);
    
        % if recessive lethal drive
        if recessivelethal
            new_larva(:,10,1)=0;
        end
    
        if haplolethal
            new_larva(:,[4,7,9],1)=0;
        end
    
        all_numbers=sum(new_larva(:,:,1),2);
        old_numbers=sum(sum(sum(mosquitos(:,:,:,2),3),4),2);
        competition=(all_numbers*num_eggs*reproduction_rate+old_larva_competition_factor*old_numbers)./(expected_larva+expected_old*old_larva_competition_factor);
        new_larva=(expected_old*num_eggs*reproduction_rate./expected_larva).*ones(n,1,1).*new_larva.*((Low_density_growth_rate./((Low_density_growth_rate-1)*competition+1)).*ones(n,1,1));
        new_larva(:,2,1)=new_larva(:,2,1)*drive_hete_fitness;
        new_larva(:,5,1)=new_larva(:,5,1)*drive_fitness;
        new_larva(:,6,1)=new_larva(:,6,1)*drive_hete_fitness;
        new_larva(:,7,1)=new_larva(:,7,1)*drive_hete_fitness;
    
        mosquitos(:,:,:,1)=new_larva;
    
        %% 
        %{
        % new_larva=zeros(n,10,3);
        if (i>500-decrease_speed&&i<500+decrease_speed)
            new_larva(:,1,1)=new_larva(:,1,1)+(1-0.9*(i-500+decrease_speed)/decrease_speed/2);
        elseif(i>=500&&i<1000)
            new_larva(:,1,1)=new_larva(:,1,1)+0.1;
        else
            new_larva(:,1,1)=new_larva(:,1,1)+1;
        end
        mosquitos(:,:,:,1)=new_larva;
        %}
        %get malaria
        mosquito_infectious=sum(sum(mosquitos(:,:,3,3:8),4),2)-sum(sum(mosquitos(:,5,3,3:8),4),2)*reducemth1-sum(sum(mosquitos(:,[2,6,7],3,3:8),4),2)*reducemth2;
        
        %infection
        if findmosquito
            prevalence=prevalences;
        else
            prevalence=humans(:,4);
        end

        
        mosquitos(:,:,2,3:8)=mosquitos(:,:,2,3:8)+mosquitos(:,:,1,3:8).*prevalence.*htom;
        mosquitos(:,[2,6,7],2,3:8)=mosquitos(:,[2,6,7],2,3:8)-mosquitos(:,[2,6,7],1,3:8).*prevalence.*htom*reducehtm1;
        mosquitos(:,5,2,3:8)=mosquitos(:,5,2,3:8)-mosquitos(:,5,1,3:8).*prevalence.*htom*reducehtm2;
        
        mosquitos(:,:,1,3:8)=mosquitos(:,:,1,3:8)-mosquitos(:,:,1,3:8).*(prevalence.*htom*reproduction_rate);
        mosquitos(:,[2,6,7],1,3:8)=mosquitos(:,[2,6,7],1,3:8)+mosquitos(:,[2,6,7],1,3:8).*(prevalence.*htom*reproduction_rate*reducehtm1);
        mosquitos(:,5,1,3:8)=mosquitos(:,5,1,3:8)+mosquitos(:,5,1,3:8).*(prevalence.*htom*reproduction_rate*reducehtm2);
        
        %develop
        mosquitos(:,[1,3,4,8,9,10],3,3:8)=mosquitos(:,[1,3,4,8,9,10],2,3:8)*develop_rate+mosquitos(:,[1,3,4,8,9,10],3,3:8);
        mosquitos(:,[1,3,4,8,9,10],2,3:8)=mosquitos(:,[1,3,4,8,9,10],2,3:8).*(1-develop_rate);
        mosquitos(:,5,3,3:8)=mosquitos(:,5,2,3:8)*develop_rate*(1-reducedevelop1)+mosquitos(:,5,3,3:8);
        mosquitos(:,5,2,3:8)=mosquitos(:,5,2,3:8).*(1-develop_rate.*(1-reducedevelop1));
        mosquitos(:,[2,6,7],3,3:8)=mosquitos(:,[2,6,7],2,3:8)*develop_rate*(1-reducedevelop2)+mosquitos(:,[2,6,7],3,3:8);
        mosquitos(:,[2,6,7],2,3:8)=mosquitos(:,[2,6,7],2,3:8).*(1-develop_rate.*(1-reducedevelop2));
        
        %human malaria
        %old infection way ,may cause problem
        %humans(:,2)=humans(:,1).*mosquito_infectious.*(1-human_resistance).*humanlist.mtoh;
        %humans(:,1)=1-humans(:,2)-humans(:,3)-humans(:,4);
        if findhuman
            humans(:,1)=exp(-(EIRS./52).*(humanlist.b1+(1-humanlist.b1)./(1+human_resistance.^humanlist.shape)).*humanlist.mtoh+log(humans(:,1)))+humans(:,4)*humanlist.human_recovery;
        else
            humans(:,1)=exp(-mosquito_infectious*reproduction_rate.*(humanlist.b1+(1-humanlist.b1)./(1+human_resistance.^humanlist.shape)).*humanlist.mtoh+log(humans(:,1)))+humans(:,4)*humanlist.human_recovery;
        end
        humans(:,4)=humans(:,4)*(1-humanlist.human_recovery)+humans(:,3);
        humans(:,3)=humans(:,2);
        humans(:,2)=1-humans(:,1)-humans(:,3)-humans(:,4);
    
        %prevent very low density human transmission
        %humans(humans(:,[1,4])<0.01)=0;
        %humans(:,4)=0;
        
        %human_resistance=human_resistance+humans(:,4)*humanlist.immunity_gain_rate-(1-humans(:,4))*humanlist.immunity_losing_rate;
        human_resistance=human_resistance+humans(:,4)*humanlist.immunity_gain_rate-human_resistance*humanlist.immunity_losing_rate;
    
    
        %store output
        %{
        mosquito_malarias(:,i)=mosquito_infectious;
        human_malarias(:,i)=humans(:,4);
        human_immunitys(:,i)=human_resistance;
        %}
    
        %plot output
        gene2=squeeze(2*genotypes(:,5)+genotypes(:,2)+genotypes(:,6)+genotypes(:,7));
        gene1=squeeze(2*genotypes(:,1)+genotypes(:,2)+genotypes(:,3)+genotypes(:,4));
        gene3=squeeze(2*genotypes(:,8)+genotypes(:,6)+genotypes(:,3)+genotypes(:,9));
        gene4=squeeze(2*genotypes(:,10)+genotypes(:,9)+genotypes(:,4)+genotypes(:,7));
        if plotgif
            clf(f1);
    
        
            subplot(2,1,1),    
            plot(x,max,'black', x, gene1, 'blue', x, gene2,x, gene3, x, gene4);
            title('mosquito gene frequency')
            legend('','wild type','drive','functional resistance','nonfunctional resistance');
    
            subplot(2,1,2)

            if findmosquito
                plot(mosquito_infectious*reproduction_rate*52,prevalences,'green');
            
            else
                plot( x, humans(:,4), 'red' , x,1-(humanlist.b1+(1-humanlist.b1)./(1+human_resistance.^humanlist.shape)),'blue', x,mosquito_infectious*reproduction_rate/7,'green');
            end
            result = min(mean(humans(:,4)), result);
            EIRS1=[3.748022899,4.209665958,20.302014,19.46226874,37.85310084,68.37694864,41.62689645,67.65891383,92.87530606,128.8428103,126.1500306,157.4671712,198.6464579,209.4149771,240.2292953,169.5482093,159.1395626,162.5352447];
            EIRS1=reshape(EIRS1,length(EIRS1),1);
            prevalences1=[0.02769617,0.167517702,0.098683032,0.283678023,0.25356346,0.268620741,0.612798931,0.511697506,0.470825206,0.685937987,0.617101703,0.59559107,0.580532176,0.552568192,0.602044422,0.806399466,0.892445223,0.944072032];
            prevalences1=reshape(prevalences1,length(prevalences1),1);
            hold on
            % scatter(EIRS1,prevalences1,'.');
            title('malaria')
            legend('human prelalance','immunity','daily EIR');%,'vector_density'
        
        
            make_gif(storename, i);
        end
        
        %fussion
        if n>=e+1
            new_mosquito=mosquitos(:,:,:,3:8);
            for p=1:e
                new_mosquito(p+1:n,:,:,:)=new_mosquito(p+1:n,:,:,:)+mosquitos(1:n-p,:,:,3:8).*fussion(p)-mosquitos(p+1:n,:,:,3:8).*fussion(p);
                new_mosquito(1:n-p,:,:,:)=new_mosquito(1:n-p,:,:,:)+mosquitos(p+1:n,:,:,3:8).*fussion(p)-mosquitos(1:n-p,:,:,3:8).*fussion(p);
            end
            mosquitos(:,:,:,3:8)=new_mosquito;
        end

        if sample
            if gene1(samplepoint1)<capacity
                if i<t1
                    t1=i;
                end
            end
            if gene1(samplepoint2)<capacity
                t2=i;
                break
            end
        end
    
        
    end
    if sample
        result= (samplepoint2-samplepoint1)/(t2-t1)/speed;
        
    end
    if findhuman
        delta=prevalences-humans(:,4);
        result=sum(delta.*delta);
    end
    if findmosquito
        %mosquito_infectious*reproduction_rate*capacity*52
        delta=EIRS-reproduction_rate.*mosquito_infectious.*52;
        result=sum(delta.*delta);
    end

    

