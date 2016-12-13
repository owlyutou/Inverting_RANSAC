function Hout = robustfit_Homog(Hin,x1,x2,nInliers,wFunName)
if ~exist('wFunName','var')
    wFunName = 'bisquare';
end
assert(all(x1(3,:)==1))
assert(all(x2(3,:)==1))
if size(Hin,1) == 9
    Hin = reshape(Hin,3,3)';
end

%% sort by error and take inliers
trans = @(h,x)bsxfun(@rdivide,h*x,h(3,:)*x);
err   = sum((trans(Hin,x1)-x2).^2);
[~,perm] = sort(err);

x1_in = x1(:,perm(1:nInliers));
x2_in = x2(:,perm(1:nInliers));

%% calc robust H

vec_1 = ones( nInliers,1);
vec_0 = zeros(nInliers,3);

U = x2_in(1:2,:)';

X = [x1_in(1:2,:)'   vec_1  vec_0                  [-x2_in(1,:).*x1_in(1,:) ; -x2_in(1,:).*x1_in(2,:)]';
    vec_0                   x1_in(1:2,:)'   vec_1  [-x2_in(2,:).*x1_in(1,:) ; -x2_in(2,:).*x1_in(2,:)]'  ];

Tvec = zeros(9,1);
% Tvec(1:8) = X \ U;

% this is a slightly different version of Matlab's "robustfit", where
% samples are rewighted in pairs (xy coordinates). The original version 
% will also work well (probably) but is not "correct". (Roee)
Tvec(1:8) = robustfit_local(X , U(:),wFunName,[],0);


% We assumed I = 1;
Tvec(9) = 1;

Hout = Tvec; % reshape(Tvec,3,3)';
