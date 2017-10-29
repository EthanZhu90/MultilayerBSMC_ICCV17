function Vxx = bsmc_computeVxx(w, h)

n = w * h;

I = []; J = []; V=[];
% Connect vertically (down, then up)
is = (1:n)'; is(h:h:n)=[];
js = is+1;
I = [I;is;js];
J = [J;js;is];
V = [V;-1*ones(size(is)); -1*ones(size(js))];

% Connect horizontally (right, then left)
is = (1:(n-h))';
js = is+h;
I = [I;is;js];
J = [J;js;is];
V = [V;-1*ones(size(is)); -1*ones(size(js))];

% Add diagonal wieghts
% For all nonboundary points they are connected to 4 pixels 
[x,y] = meshgrid(2:(w-1), 2:(h-1));
is = sub2ind([h w], y(:), x(:));
I = [I;is];
J = [J;is];
V = [V; 4*ones(size(is))];
% For nonboundary pixels they are connected to 3 pixels each except the 4
% corners they are connected to 2.

% Add three connected pixels (diagonal entries)
is = [2:(h-1) ...
    (h+1):h:((h*(w-2))+1) ...
    2*h:h:h*(w-1) ...
    (h*(w-1)+2):((h*w)-1)]';
I = [I; is];
J = [J; is];
V = [V; 3*ones(size(is))];

% Add two connected pixels 
is = [1 h (h * (w-1))+1 h*w]';
I = [I;is];
J = [J;is];
V = [V; 2*ones(size(is))];
Vxx = sparse(I,J,V,n,n);
end