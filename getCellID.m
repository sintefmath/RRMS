function velocity = getVelocity(G,states, p)
id = getCellID(G,p)
velocity = 
end
function id = getCellID(G,p)
id=NaN;
for i=1:G.cells.num
    
    [n, pn]=gridCellNodes(G,i);
    %assuming tetrahedral grid from now on
    a=G.nodes.coords(n(1),:);
    b=G.nodes.coords(n(2),:);
    c=G.nodes.coords(n(3),:);
    d=G.nodes.coords(n(4),:);
    tmp=isinside(a,b,c,d, p);
    if tmp
        disp(strcat(num2str(i), 'yes'))
    %else
    %    disp(strcat(num2str(i), 'no'))
    end
    if tmp
        id=i;
        return
    end
end



function val = isinside(a,b,c,d, p)

vap = p - a;
vbp = p - b;

vab = b - a;
vac = c - a;
vad = d - a;

vbc = c - b;
vbd = d - b;

va6 = dot(cross(vbp, vbd), vbc);
vb6 = dot(cross(vap, vac), vad);
vc6 = dot(cross(vap, vad), vab);
vd6 = dot(cross(vap, vab), vac);
v6 = 1 / dot(cross(vab, vac), vad);
bary = [va6*v6, vb6*v6, vc6*v6, vd6*v6];
if max(bary)<=1 && min(bary)>=0
    val = true;
else
    val = false;
end