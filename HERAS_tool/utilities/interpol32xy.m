function Plinexy = interpol32xy(F,uv,V,Pline)
% this function takes a map (tri, uv and xy) and a point in uv, and
% calculates its position in xy
% remove the nan
% Authors : Salvador Mercader Pellicer
isNAN = find(isnan(Pline(:,1)));
Pline(isNAN,:) = zeros(size(Pline(isNAN,:)));

[t,P] = tsearchn(uv',F',Pline);
Plinexy = nan(size(P,1),3);
for ii = 1:length(t)
    if ~isnan(t(ii))
        PPtemp = V(:,F(:,t(ii))).*(P(ii,:).*ones(3,1));
        Plinexy(ii,:) = sum(PPtemp,2);
    end
end

Plinexy(isNAN,:) = nan(size(Plinexy(isNAN,:)));
end
