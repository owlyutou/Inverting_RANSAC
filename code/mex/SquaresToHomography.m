function [H] = SquaresToHomography(sq1x,sq1y,sq2x,sq2y)
% sq1x, sq1y are double: 4x1
% sq2x, sq2y are int16: 4xN

if coder.target('MATLAB')
    % Code for MATLAB execution, ignored by Matlab coder
    if exist('SquaresToHomography_mex','file')==3
        H = SquaresToHomography_mex(sq1x,sq1y,sq2x,sq2y);
        return
    else
        warning('SquaresToHomography:NoMex','!!! no mex file found, running matlab version is MUCH slower!!!')
    end
end


[tmp1,N1] = size(sq1x);
[tmp2,N2] = size(sq2x);

assert(4== tmp1 && tmp1 == tmp2)
assert(N1==1)

assert(tmp1 == size(sq1y,1) && N1 == size(sq1y,2))
assert(tmp2 == size(sq2y,1) && N2 == size(sq2y,2))

H = zeros(9,N2);
for ii = 1:N2
    H(1:9,ii) = CreatreOneHomography(sq1x,sq1y,sq2x(1:4,ii),sq2y(1:4,ii));
end



function H = CreatreOneHomography(sq1x,sq1y,sq2x,sq2y)


vec_1 = ones( 4,1);
vec_0 = zeros(4,3);


U = double([sq2x; sq2y]);

X = [sq1x   sq1y   vec_1  vec_0                -double(sq2x).*sq1x  -double(sq2x).*sq1y;
     vec_0                sq1x   sq1y   vec_1  -double(sq2y).*sq1x  -double(sq2y).*sq1y  ];

Tvec = zeros(9,1);
Tvec(1:8) = X \ U;

% We assumed I = 1;
Tvec(9) = 1;

H = Tvec;
% H = reshape(Tvec,3,3);

