function matRad_sliderUI(data,dij,pln,ct)
    fInds = data.finds;
    
    
    %sort function values according to each dimension
    [A,I] = sort(fInds(:,1),1,'ascend');
    %%
    fIndsSorted = fInds(I,:);
    
    %%
    
    weights = data.weights(:,I);
    %%
    pics = {};
    
    slice = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);
    
    for i = 1:size(weights,2)
        resultGUI = matRad_calcCubes(weights(:,i),dij);
        pics{i} = resultGUI.physicalDose(:,:,slice);
    end
    %%
    sortIdxs = zeros(size(fInds));
    sliderIdxs = zeros(size(fInds)); %stores the slider combinations that belong together
    sliderMapper = zeros(size(fInds)); 
    %%
    for j = 1:size(fIndsSorted,2)
        [A,I] = sort(fIndsSorted(:,j),1,'ascend');
        sortIdxs(:,j) = I;
    end
    %%
    for i = 1:size(fIndsSorted,1)
        [row,col] = find(sortIdxs==i);
        sliderIdxs(i,:) = row;
    end
    %%
    %%
    f = figure;
    
    axes('Position',[.4 .5 .4 .4]);
    n = size(fInds,1);
    p = uipanel(f,'Position',[0.1 0.1 0.20 0.35]);
    
    names = {'Objective 3','Objective 2','Objective 1'};
    slidInit = 1;
    sliders = {};
    namesui = {};
    sliderValues = {};
    
    for i = 1:size(fInds,2)
        namesui{i} = uicontrol(p,'Style','text',...
            'Position',[5 10+(i-1)*25 70 20],...
            'String',names{end-i+1});
    
        sliderValues{i} = uicontrol(p,'Style','text',...
            'Position',[220 10+(i-1)*25 70 20],...
            'String',fIndsSorted(1,i));
    
        sliders{i} = uicontrol(p,'Style','slider',...
            'Min',1,'Max',n,...
            'SliderStep', [1/(n-1) 1],...
            'Position',[100 10+(i-1)*25 100 20]);
                
        sliders{i}.Value = sortIdxs(1,i);
    
    end 
    
    for i = 1:size(fInds,2)
        set(sliders{i},'Callback',{@slider_callback,sliderValues,sliders,fIndsSorted,sortIdxs,sliderIdxs,pics,i})
    end
end
    % End of main file
    
    % Callback subfunctions to support UI actions
    function slider_callback(slider,~,textFields,sliders,fIndsSorted,sortIdxs,sliderIdxs,pics,idx)
        % Update the text for the moved slider
        set(textFields{idx},'String',fIndsSorted(sortIdxs(slider.Value,idx),idx))
    
    
        for i = setdiff(1:numel(sliders),idx)
            set(sliders{i},'Value',sliderIdxs(sortIdxs(slider.Value,idx),i))
            set(textFields{i},'String',fIndsSorted(sortIdxs(slider.Value,idx),i))
        end
        plotterfcn(pics{sortIdxs(slider.Value,idx)})
    end
    
    function plotterfcn(pic)
        % Plots the image
        imagesc(pic),colorbar, colormap(jet);
    end
