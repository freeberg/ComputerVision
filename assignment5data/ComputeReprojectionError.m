function [err,res] = ComputeReprojectionError(P,U,u)
%Compute the reprojection error from the current solution P,U and the
%imagedata u. The value of each residual is in res.

err = 0;
res = [];
for i = 1:length(P);
    uu = u{i};
    vis = isfinite(uu(1,:));
    err = err + ...
        sum(((P{i}(1,:)*U(:,vis))./(P{i}(3,:)*U(:,vis)) - uu(1,vis)).^2) + ...
        sum(((P{i}(2,:)*U(:,vis))./(P{i}(3,:)*U(:,vis)) - uu(2,vis)).^2);
    res = [res ((P{i}(1,:)*U(:,vis))./(P{i}(3,:)*U(:,vis)) - uu(1,vis)).^2 + ...
            ((P{i}(2,:)*U(:,vis))./(P{i}(3,:)*U(:,vis)) - uu(2,vis)).^2];
end
