%Grafica los parámetros en el tiempo y ordenado por correlación. Cada punto
%representa un hamster y cada franja una marca temporal.

clc, clear all ,close all;

%% Parámetros Iniciales.
Grupo='C';%Grupo de Longitudes
numFechas=1:12;%Rango de Fechas
T=1;%Tipo de Leishmaniasis
E=1;%Estado - Criterio de curación
linux=0;
compararFirmaControl=0;

%% Parámetros Constantes
if compararFirmaControl
    load('meanControlHSS.mat');
end
tipos={'LB-TE','LB-TC','LB-ST','LP-TE','LP-TC','LP-ST'};
estados={'CURADO','MEJORIA','NOCURADO'};
areas={'SANO','BORDE','CENTRO'};
VariablesNames={'Epidermis','EpiPDermis','Keratynocites','Collagen',...
    'FColl','FibroblastD','MacrophagesD','Fmel','FBlood','OS'};
colorName=[{'Rojo'},{'Verde'},{'Azul'},{'Negro'},{'Cyan'}];
load('Datos');
dataB=struct;
dataC=struct;
%% Carga de Datos
Dir='/media/labmirp/OS/Users/Labmirp/Documents/RicardoF/Proyecto Investigacion/FirmasActualizadas/';
% Datos por Grupo
if linux
    datosG1=load(strcat('/media/labmirp/OS/Users/Labmirp/Documents/RicardoF/Proyecto Investigacion/FirmasActualizadas/G1',Grupo,'/Total.mat'));
    datosG2=load(strcat('/media/labmirp/OS/Users/Labmirp/Documents/RicardoF/Proyecto Investigacion/FirmasActualizadas/G2',Grupo,'/Total.mat'));
    datosG3=load(strcat('/media/labmirp/OS/Users/Labmirp/Documents/RicardoF/Proyecto Investigacion/FirmasActualizadas/G3',Grupo,'/Total.mat'));
    datosG4=load(strcat('/media/labmirp/OS/Users/Labmirp/Documents/RicardoF/Proyecto Investigacion/FirmasActualizadas/G4',Grupo,'/Total.mat'));
else
    datosG1=load(strcat('C:\Users\Labmirp\Documents\RicardoF\Proyecto Investigacion\FirmasActualizadas\G1',Grupo,'\Total.mat'));
    datosG2=load(strcat('C:\Users\Labmirp\Documents\RicardoF\Proyecto Investigacion\FirmasActualizadas\G2',Grupo,'\Total.mat'));
    datosG3=load(strcat('C:\Users\Labmirp\Documents\RicardoF\Proyecto Investigacion\FirmasActualizadas\G3',Grupo,'\Total.mat'));
    datosG4=load(strcat('C:\Users\Labmirp\Documents\RicardoF\Proyecto Investigacion\FirmasActualizadas\G4',Grupo,'\Total.mat'));
end

%% Proceso

% Crea estrucutras ordenadas por fecha y correlación
for Date=numFechas % Recorre Fechas
    %disp(Date)
    dataSortB={};
    dataSortC={};
    dataSano=zeros(1961,1);
    meanBorde=inf*ones(1961,1);
    meanCentro=inf*ones(1961,1);
    for t=tipos(T)
        data=Datos(strcmp({Datos.Tipo},char(t))); % Identifica los hasmter de un tipo específico
        
        for e=estados(E)
            dataEstado=data(strcmp({data.Estado},char(e)) & [data.Enabled]); % Identifica los hamster de un estado específico
            dataEstado=assingColor(dataEstado); % Asigna un color a cada hamster
            if ~isempty(dataEstado)
                %% Recorre dataEstados
                for dE=dataEstado
                    % Carga los valores de cada Hamster
                    switch dE.Grupo
                        case 'G1'
                            valores=datosG1.TOTAL(strcmp({datosG1.TOTAL.NAME},dE.Hamster));
                        case 'G2'
                            valores=datosG2.TOTAL(strcmp({datosG2.TOTAL.NAME},dE.Hamster));
                        case 'G3'
                            valores=datosG3.TOTAL(strcmp({datosG3.TOTAL.NAME},dE.Hamster));
                        case 'G4'
                            valores=datosG4.TOTAL(strcmp({datosG4.TOTAL.NAME},dE.Hamster));
                    end
                    
                    % Elimina fechas vacias
                    tab=cell2table({valores.VART.Date}');
                    idx=all(cellfun(@isempty,tab{:,:}),2);
                    valores.VART=valores.VART(~idx);
                    
                    %Si hay valores de Borde y Centro halla la firma mediana, si no asume correlación de infinito.
                    try
                        val= valores.VART(Date); %Toma valores específicos de la fecha
                        optSano=sum(val.SANO_OK==1);
                        dataSano= [dataSano,reshape(val.SANO_var(:,2,val.SANO_OK),[1961,optSano])];
                        try
                            optBorde=sum(val.BORDE_OK==1);
                            meanBorde=[meanBorde,mean(reshape(val.BORDE_var(:,2,val.BORDE_OK),[1961,optBorde]),2)];
                        catch
                            meanBorde=[meanBorde,inf*ones(1961,1)];
                        end
                        try
                            optCentro=sum(val.CENTRO_OK==1);
                            meanCentro=[meanCentro,mean(reshape(val.CENTRO_var(:,2,val.CENTRO_OK),[1961,optCentro]),2)];
                        catch
                            meanCentro=[meanCentro,inf*ones(1961,1)];
                        end
                    catch
                    end
                    
                end
                %% Corregir Datos
                
                dataSano=dataSano(:,2:end);
                meanCentro=meanCentro(:,2:end);
                meanBorde=meanBorde(:,2:end);
                meanSano=mean(dataSano,2);
                
                
                %% Correlacion - Centro
                % Primera Fila
                con2=0;
                con=1;
                try
                    if compararFirmaControl
                        corr=corrcoef(meanHSS,meanCentro(:,con));
                    else
                        corr=corrcoef(meanSano,meanCentro(:,con));
                    end
                    
                    
                    if ~isnan(sum(corrcoef(meanSano,meanCentro(:,con))))
                        con2=con2+1;
                        coefcorrC=min(min(corr));
                        dataSortC={dataEstado(1).Tipo,dataEstado(1).Estado,dataEstado(1).Grupo,dataEstado(1).Hamster, coefcorrC,dataEstado(1).Color,dataEstado(1).ColorName};
                    end
                catch
                end
                %Filas 2 a fin
                
                for dE=dataEstado(2:end)
                    con=con+1;
                    try
                        if compararFirmaControl
                            corr=corrcoef(meanHSS,meanCentro(:,con));
                        else
                            corr=corrcoef(meanSano,meanCentro(:,con));
                        end
                        
                        if ~isnan(sum(corrcoef(meanSano,meanCentro(:,con))))
                            con2=con2+1;
                            coefcorrC=min(min(corr));
                            dataSortC{con2,1}=dE.Tipo;
                            dataSortC{con2,2}=dE.Estado;
                            dataSortC{con2,3}=dE.Grupo;
                            dataSortC{con2,4}=dE.Hamster;
                            dataSortC{con2,5}=coefcorrC;
                            dataSortC{con2,6}=dE.Color;
                            dataSortC{con2,7}=dE.ColorName;
                        end
                    catch
                    end
                end
                
                %% Ordenar Datos - Centro
                try
                    dataSortC=sortrows(dataSortC,[-5,1]);
                    con=1;
                    for ds=dataSortC'
                        dataC(end+1).Grupo=ds(3);
                        dataC(end).Hamster=ds(4);
                        dataC(end).Date=Date;
                        dataC(end).Pos=con;
                        dataC(end).Tipo=t;
                        dataC(end).Estado=e;
                        dataC(end).Correlacion=ds{5};
                        dataC(end).Color=ds{6};
                        dataC(end).ColorName=ds{7};
                        con=con+1;
                    end
                catch
                end
                
                %% Correlacion - Borde
                % Primera Fila
                con2=0;
                con=1;
                try
                    if compararFirmaControl
                        corr=corrcoef(meanHSS,meanBorde(:,con));
                    else
                        corr=corrcoef(meanSano,meanBorde(:,con));
                    end
                    if ~isnan(sum(corrcoef(meanSano,meanBorde(:,con))))
                        con2=con2+1;
                        coefcorrB=min(min(corr));
                        dataSortB={dataEstado(1).Tipo,dataEstado(1).Estado,dataEstado(1).Grupo,dataEstado(1).Hamster, coefcorrB,dataEstado(1).Color,dataEstado(1).ColorName};
                    end
                catch
                end
                %Filas 2 a fin
                
                for dE=dataEstado(2:end)
                    con=con+1;
                    try
                        if compararFirmaControl
                            corr=corrcoef(meanHSS,meanBorde(:,con));
                        else
                            corr=corrcoef(meanSano,meanBorde(:,con));
                        end
                        if ~isnan(sum(corrcoef(meanSano,meanBorde(:,con))))
                            con2=con2+1;
                            coefcorrB=min(min(corr));
                            dataSortB{con2,1}=dE.Tipo;
                            dataSortB{con2,2}=dE.Estado;
                            dataSortB{con2,3}=dE.Grupo;
                            dataSortB{con2,4}=dE.Hamster;
                            dataSortB{con2,5}=coefcorrB;
                            dataSortB{con2,6}=dE.Color;
                            dataSortB{con2,7}=dE.ColorName;
                        end
                    catch
                    end
                end
                
                %% Ordenar Datos - Borde
                try
                    dataSortB=sortrows(dataSortB,[-5,1]);
                    con=1;
                    for ds=dataSortB'
                        dataB(end+1).Grupo=ds(3);
                        dataB(end).Hamster=ds(4);
                        dataB(end).Date=Date;
                        dataB(end).Pos=con;
                        dataB(end).Tipo=t;
                        dataB(end).Estado=e;
                        dataB(end).Correlacion=ds{5};
                        dataB(end).Color=ds{6};
                        dataB(end).ColorName=ds{7};
                        con=con+1;
                    end
                catch
                end
                
                
                
            end
        end
    end
end

%% Recorre las Fechas y Grafica los Datos

% Corrige los datos
tableB=struct2table(dataB(2:end))
tableC=struct2table(dataC(2:end))

fC=figure;
fB=figure;
save('TablesData.mat','tableC','tableB');
for t=tipos(T)
    data=Datos(strcmp({Datos.Tipo},char(t)));
    for e=estados(E)
        dataEstado=data(strcmp({data.Estado},char(e)) & [data.Enabled]);
        dataEstado=assingColor(dataEstado); % Asigna un color a cada hamster
        
        %dataEstado=data(strcmp({data.Estado},char(e)));
        con=1;
        for dE = dataEstado
            name=strcat(dE.Grupo,'_',dE.Hamster);
            try
                tE=tableB(strcmp(tableB.Hamster,dE.Hamster) & strcmp(tableB.Grupo,dE.Grupo),:);
                figure(fB);
                
                plot(tE.Date,tE.Correlacion,'Color',[dE.Color])
                title(strcat('Correlaciones Borde Firmas',Grupo));
                xlabel('Fechas');
                ylabel('Correlaci�n')
                grid on;
                hold on;
                scatter(tE.Date,tE.Correlacion,[],[dE.Color],'filled')
                hold on;
            catch
            end
            tE=tableC(strcmp(tableC.Hamster,dE.Hamster) & strcmp(tableC.Grupo,dE.Grupo),:);
            figure(fC)
            
            plot(tE.Date,tE.Correlacion,'Color',[dE.Color])
            title(strcat('Correlaciones Centro Firmas',Grupo));
            xlabel('Fechas');
            ylabel('Correlaci�n')
            grid on;
            hold on;
            scatter(tE.Date,tE.Correlacion,[],[dE.Color],'filled')
            hold on;
            
            con=con+1;
        end
    end
end
%% GRaf PARAMs
figure;
for t=tipos(T)
    for e=estados(E)
        for var=1:size(VariablesNames,2)
            
            %% Grafica Datos - Borde
            subplot(2,1,1),
            data=Datos(strcmp({Datos.Tipo},char(t)));
            dataEstado=data(strcmp({data.Estado},char(e)) & [data.Enabled]);
            dataEstado=assingColor(dataEstado);
            y=[];
            datosPorFecha=size(dataEstado,2);
            newXLabel=0:0.25:1;
            currentDataB= tableB(strcmp(tableB.Tipo,t),:);
            currentDataB= currentDataB(strcmp(currentDataB.Estado,e),:);
            line=struct;
            for date=numFechas
                dateData=currentDataB(currentDataB.Date==date,:);
                if(~isempty(dateData))
                    s=size(dateData,1);
                    for j=1:s
                        dt=dateData(j,:);
                        switch char(dt.Grupo)
                            case 'G1'
                                valores=datosG1.TOTAL(strcmp({datosG1.TOTAL.NAME},char(dt.Hamster)));
                            case 'G2'
                                valores=datosG2.TOTAL(strcmp({datosG2.TOTAL.NAME},char(dt.Hamster)));
                            case 'G3'
                                valores=datosG3.TOTAL(strcmp({datosG3.TOTAL.NAME},char(dt.Hamster)));
                            case 'G4'
                                valores=datosG4.TOTAL(strcmp({datosG4.TOTAL.NAME},char(dt.Hamster)));
                        end
                        try
                            v=valores.VART(date);
                            if(~(isempty(v.SANO_AG)))
                                if(sum(v.SANO_OK)==1)
                                    medianS=v.SANO_AG(v.SANO_OK,:);
                                else
                                    medianS=median(v.SANO_AG(v.SANO_OK,:));
                                end
                                
                                y=[y,medianS(var)];
                                posx=(dt.Date-1)+(dt.Date-1)*0.25+dt.Correlacion;
                                scatter(posx,medianS(var),[],dt.Color,'s')
                                hold on;
                            end
                            if(~(isempty(v.BORDE_AG)))
                                if(sum(v.BORDE_OK)==1)
                                    medianB=v.BORDE_AG(v.BORDE_OK,:);
                                else
                                    medianB=median(v.BORDE_AG(v.BORDE_OK,:));
                                end
                                
                                y=[y,medianB(var)];
                                posx=(dt.Date-1)+(dt.Date-1)*0.25+dt.Correlacion;
                                line(end+1).posX=posx;
                                line(end).median=medianB(var);
                                line(end).colorName=dt.ColorName;
                                line(end).color=dt.Color;
                                
                                
                                scatter(posx,medianB(var),[],dt.Color,'filled')
                                hold on;
                            end
                        catch
                        end
                    end
                    
                end
            end
            line=line(2:end);
            if ~isempty(line)
                for colors=colorName(1:datosPorFecha)
                    grafData=line(strcmp([line.colorName],colors)');
                    plot([grafData.posX],[grafData.median],'Color',grafData(1).color)
                end
            end
            set(gca,'xtick',0:0.25:15);
            set(gca,'XtickLabel',newXLabel);
            tam=max(ylim)*1.2;
            for i=min(numFechas):max(numFechas)
                if mod(i,2)
                    m2=(i)+0.125*((i*2)-1);
                    m1=m2-1.25;
                    box1=[m1 m1  m2 m2];
                    boxy=[0 tam tam 0];
                    patch(box1,boxy,[0 0 1],'FaceAlpha',0.2)
                end
            end
            grid on;
            try
                ylim(1.1*[0 max(y)]);
                x=xlim;
                xlim([-0.125 ,x(2)]);
            catch
            end
            title(strcat('BORDE - ',VariablesNames(var),'-',t,'-',e));
            
            %% Grafica Datos - Centro
            subplot(2,1,2),
            
            data=Datos(strcmp({Datos.Tipo},char(t)));
            dataEstado=data(strcmp({data.Estado},char(e)) & [data.Enabled]);
            dataEstado=assingColor(dataEstado);
            y=[];
            datosPorFecha=size(dataEstado,2);
            newXLabel=0:0.25:1;
            currentDataC= tableC(strcmp(tableC.Tipo,t),:);
            currentDataC= currentDataC(strcmp(currentDataC.Estado,e),:);
            line=struct;
            for date=numFechas
                dateData=currentDataC(currentDataC.Date==date,:);
                if(~isempty(dateData))
                    s=size(dateData,1);
                    for j=1:s
                        dt=dateData(j,:);
                        switch char(dt.Grupo)
                            case 'G1'
                                valores=datosG1.TOTAL(strcmp({datosG1.TOTAL.NAME},char(dt.Hamster)));
                            case 'G2'
                                valores=datosG2.TOTAL(strcmp({datosG2.TOTAL.NAME},char(dt.Hamster)));
                            case 'G3'
                                valores=datosG3.TOTAL(strcmp({datosG3.TOTAL.NAME},char(dt.Hamster)));
                            case 'G4'
                                valores=datosG4.TOTAL(strcmp({datosG4.TOTAL.NAME},char(dt.Hamster)));
                        end
                        %try
                        v=valores.VART(date);
                        if(~(isempty(v.SANO_AG)))
                            if(sum(v.SANO_OK)==1)
                                medianS=v.SANO_AG(v.SANO_OK,:);
                            else
                                medianS=median(v.SANO_AG(v.SANO_OK,:));
                            end
                            y=[y,medianS(var)];
                            posx=(dt.Date-1)+(dt.Date-1)*0.25+dt.Correlacion;
                            scatter(posx,medianS(var),[],dt.Color,'s')
                            hold on;
                        end
                        if(~(isempty(v.CENTRO_AG)))
                            if(sum(v.CENTRO_OK)==1)
                                medianC=v.CENTRO_AG(v.CENTRO_OK,:);
                            else
                                medianC=median(v.CENTRO_AG(v.CENTRO_OK,:));
                            end
                            
                            y=[y,medianC(var)];
                            posx=(dt.Date-1)+(dt.Date-1)*0.25+dt.Correlacion;
                            line(end+1).posX=posx;
                            line(end).median=medianC(var);
                            line(end).colorName=dt.ColorName;
                            line(end).color=dt.Color;
                            
                            scatter(posx,medianC(var),[],dt.Color,'filled')
                            hold on;
                        end
                        %catch
                        %end
                    end
                    
                end
            end
            line=line(2:end);
            if ~isempty(line)
                for colors=colorName(1:datosPorFecha)
                    grafData=line(strcmp([line.colorName],colors)');
                    plot([grafData.posX],[grafData.median],'Color',grafData(1).color)
                end
            end
            set(gca,'xtick',0:0.25:15);
            set(gca,'XtickLabel',newXLabel);
            tam=max(ylim)*1.2;
            for i=min(numFechas):max(numFechas)
                if mod(i,2)
                    m2=(i)+0.125*((i*2)-1);
                    m1=m2-1.25;
                    box1=[m1 m1  m2 m2];
                    boxy=[0 tam tam 0];
                    patch(box1,boxy,[0 0 1],'FaceAlpha',0.2)
                end
            end
            grid on;
            try
                ylim(1.1*[0 max(y)]);
                x=xlim;
                xlim([-0.125 ,x(2)]);
            catch
            end
            title(strcat('CENTRO - ',VariablesNames(var),'-',t,'-',e));
            pause;
            clf('reset')
            
        end
        
    end
    
end






