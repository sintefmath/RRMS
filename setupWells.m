%%% Copyright 2018 Equinor ASA

function [W] = setupWells(model)

mrstModule add ad-core

%% find all indices that satisfy

maxdist=10;

% vertical producers (high pressure)
w_x=0;
w_y=0;

productionInx=[];
for i=1:model.G.cells.num
    dist = (model.G.cells.centroids(i,1)-w_x).^2;
    dist = dist + (model.G.cells.centroids(i,2)-w_y).^2;
    if (sqrt(dist)<maxdist)
        productionInx(end+1) = i;
    end
end
%WI=ones(length(productionInx),1);
%W = addWell([], model.G, model.rock, productionInx, 'Sign', -1, 'Dir', 'z', 'WI', WI);% , 'Radius', 0.1, 'Comp_i', [0, 1]);
W = addWell([], model.G, model.rock, productionInx, ...
             'InnerProduct', 'ip_tpf', ...
            'Type', 'bhp', 'Val', 100.0*barsa, ...
            'Radius', 0.1, 'Comp_i', [0, 1], 'Sign', -1, 'Name', 'P');


% vertical injectors (low pressure)

w_x=700;
w_y=400;

injectionInx=[];
for i=1:model.G.cells.num
    dist = (model.G.cells.centroids(i,1)-w_x).^2;
    dist = dist + (model.G.cells.centroids(i,2)-w_y).^2;
    if (sqrt(dist)<maxdist)
        injectionInx(end+1) = i;
    end
end

WI=ones(length(injectionInx),1);
%W = addWell(W, model.G, model.rock, injectionInx, 'Sign', 1, 'Dir', 'z', 'WI', WI);% , 'Radius', 0.1, 'Comp_i', [1, 0]);
W = addWell(W, model.G, model.rock, injectionInx, ...
             'InnerProduct', 'ip_tpf', ...
            'Type', 'rate' , 'Val', 2000*meter^3/day, ...
            'Radius', 0.1, 'Dir', 'z', 'Comp_i', [1, 0], 'Sign', 1, 'Name', 'I');
