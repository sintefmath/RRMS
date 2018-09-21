function bin = fluxdistribution(tof, flux, binedges)
    assert(length(tof)==length(flux));
    bin = struct;
    numbins=length(binedges)-1;
    bin.n = zeros(numbins,1);
    bin.f = zeros(numbins,1);
    for i=1:length(tof)
        binindex = find(binedges<tof(i));
        if isempty(binindex)
            binindex = 1;
        else
            binindex = min(binindex(end), numbins);
        end
        bin.f(binindex) = bin.f(binindex)+flux(i);
    end  
    bin.f(:)./bin.n(:);
end