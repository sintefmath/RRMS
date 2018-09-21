
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
bary = [va6*v6, vb6*v6, vc6*v6, vd6*v6]
if max(bary)<=1 && min(bary)>=0
    val = true;
else
    val = false;
end
end