function out = MDVALSE( Y, L, X, D, M)
%% MDVALSE algorithm
 
[c_M{1:D}] = ndgrid(1:M);
D_set_vec = reshape(cat(D,c_M{:}),[],D);

a_theta = @(index,omega)exp(1i*(index-1)*omega');

w     = zeros(L,1);
J     = zeros(L,L);
h     = zeros(L,1);
C     = zeros(L);
T     = 500;       
A     = zeros(prod(M),L);
mse   = zeros(T,1);
Kt    = zeros(T,1);
th    = nan(L,D);
t     = 1;
G_num = 32;
R_g   = G_num*ones(1,D);
w_grid = 2*pi*(0:(G_num-1))/G_num;
Y2 = (Y(:))'*Y(:);
res   = Y;
xro = zeros(size(Y));
for l=1:L
    M_nu = min(prod(M),500);
    yI = Y(1:M_nu);
    R  = yI*yI';
    sR = zeros(M_nu-1,1);
    for i=2:M_nu
        for k=1:i-1
            sR(i-k) = sR(i-k) + R(i,k);
        end
    end
    if l==1 
        Rh  = toeplitz([sum(diag(R));sR])/M_nu;        
        evs = sort(real(eig(Rh)),'ascend');
        nu  = mean(evs(1:floor(M_nu/4)));
        K   = floor(L/2);
        rho = K/L;
        tau = (Y2/prod(M)-nu)/(rho*L);
    end
    if D>1
        res_val = abs(fftn(reshape(res,M), R_g));
    else
        res_val = abs(fft(res, R_g));
    end
    [res_max] = max(res_val(:));
    idx = ones(D,1);
    idx_all = find(res_val == res_max);
    for i = 1:D-1
        idx(D-i+1) = floor((idx_all-1)/prod(R_g(1:D-i)))+1;
        idx_all = idx_all - (idx(D-i+1)-1)*prod(R_g(1:D-i));
    end
    idx(1) = idx_all; 
    th(l,:) = w_grid(idx)-2*pi*(w_grid(idx)>pi);
    Hess = nan(D);
    for i = 1:D
        for j = 1:D
            Hess(i,j) = 2*real((a_theta(D_set_vec,th(l,: )).*(1i*D_set_vec(:,j)))'*(Y*Y')*(a_theta(D_set_vec,th(l,: )).*(1i*D_set_vec(:,i)))...
                -(a_theta(D_set_vec,th(l,:)))'*(Y*Y')*(a_theta(D_set_vec,th(l,:)).*(D_set_vec(:,i).*D_set_vec(:,j))));
        end
    end
    Hess = Hess/(nu*prod(M));
    Hess_inv = inv(Hess);
    e = eig(Hess);
    dif = prod(e<0);
    if dif
        kappa = (Ainv(exp(0.5*diag(Hess_inv))))';
    else
        kappa = kappa*0.1; 
    end
    
    A(:,l) = prod(besseli((D_set_vec-1),repmat(kappa,prod(M),1),1)./besseli(0,repmat(kappa,prod(M),1),1),2).*a_theta(D_set_vec,th(l,: ));

    w_temp = w(1:l-1); C_temp = C(1:l-1,1:l-1);
    J(1:l-1,l) = A(:,1:l-1)'*A(:,l); J(l,1:l-1) = J(1:l-1,l)'; J(l,l) = prod(M);
    h(l) = A(:,l)'*Y;
    v = nu / ( prod(M)+ nu/tau - real(J(1:l-1,l)'*C_temp*J(1:l-1,l))/nu );
    u = v .* (h(l) - J(1:l-1,l)'*w_temp)/nu;
    w(l) = u;
    ctemp = C_temp*J(1:l-1,l)/nu;
    w(1:l-1) = w_temp - ctemp*u;
    C(1:l-1,1:l-1) = C_temp + v*(ctemp*ctemp');
    C(1:l-1,l) = -v*ctemp;  C(l,1:l-1) = C(1:l-1,l)'; C(l,l) = v;
    res = Y - A(:,1:l)*w(1:l);
    
    if l==K 
        xro    = A(:,1:l)*w(1:l);
        mse(t) = norm(X-xro)^2/norm(X)^2;
        Kt(t)  = K;
    end
end

%% Start the VALSE algorithm
cont = 1;
while cont
    t = t + 1;
    [ K, s, w, C ] = maxZ( J, h, M, nu, rho, tau );
    if K>0
        nu  = real( Y2 - 2*real(h(s)'*w(s)) + w(s)'*J(s,s)*w(s) + trace(J(s,s)*C(s,s)) )/(prod(M));
        tau = real( w(s)'*w(s)+trace(C(s,s)) )/K;
        if K<L
            rho = K/L;
        else
            rho = (L-1)/L; 
        end
    else
        rho = 1/L; 
    end
    
    inz = 1:L; inz = inz(s);
    for i = 1:K
        if K == 1
            r = Y;
            eta = 2/nu * ( r * w(inz)' );
        else
            A_i = A(:,inz([1:i-1 i+1:end]));
            r = Y - A_i*w(inz([1:i-1 i+1:end]));
            eta = 2/nu * ( r * w(inz(i))' - A_i * C(inz([1:i-1 i+1:end]),i) );
        end

        [A(:,inz(i)), th(inz(i),:)] = pntFreqEst( eta,M,D_set_vec,th(inz(i),:));
    end
    J(:,s) = A(:,:)'*A(:,s);
    J(s,:) = J(:,s)';
    J(s,s) = J(s,s) - diag(diag(J(s,s))) + prod(M)*eye(K);
    h(s)   = A(:,s)'*Y;
    
    xr     = A(:,s)*w(s);
    xr0    = A(:,s)*w(s);
    mse(t) = norm(xr-X)^2/norm(X)^2;
    Kt(t)  = K;
    if (norm(xr-xro)/norm(xro)<1e-6) || (norm(xro)==0&&norm(xr-xro)==0) || (t >= T)
        cont = 0;
        mse(t+1:end) = mse(t);
        Kt(t+1:end)  = Kt(t);
    end
    xro = xr;
end
post_c = A(:,s)*C(s,s)*A(:,s)';
post_c = (post_c+post_c')/2;
out = struct('freqs',th(s,:),'amps',w(s),'x_estimate',xr,'x_estimate0',xr0,'noise_var',nu,'iterations',t,'mse',mse,'K',Kt,'post_cov',post_c);
end

function [ K, s, w, C ] = maxZ( J, h, M, nu, rho, tau )

L = size(h,1);
cnst = log(rho/(1-rho)/tau);

K = 0; 
s = false(L,1); 
w = zeros(L,1);
C = zeros(L);
u = zeros(L,1);
v = zeros(L,1);
Delta = zeros(L,1);
if L > 1
    cont = 1;
    while cont
        if K<L-1
            v(~s) = nu ./ ( prod(M) + nu/tau - real(sum(J(s,~s).*conj(C(s,s)*J(s,~s)),1))/nu );
            u(~s) = v(~s) .* ( h(~s) - J(s,~s)'*w(s))/nu;
            Delta(~s) = log(v(~s)) + u(~s).*conj(u(~s))./v(~s) + cnst;
        else
            Delta(~s) = -1; 
        end
        if ~isempty(h(s))
            Delta(s) = -log(diag(C(s,s))) - w(s).*conj(w(s))./diag(C(s,s)) - cnst;
        end
        [~, k] = max(Delta);
        if Delta(k)>0
            if s(k)==0 
                w(k) = u(k);
                ctemp = C(s,s)*J(s,k)/nu;
                w(s) = w(s) - ctemp*u(k);
                C(s,s) = C(s,s) + v(k)*(ctemp*ctemp');
                C(s,k) = -v(k)*ctemp;
                C(k,s) = C(s,k)';
                C(k,k) = v(k);
                s(k) = ~s(k); K = K+1;
            else 
                s(k) = ~s(k); K = K-1;
                w(s) = w(s) - C(s,k)*w(k)/C(k,k);
                C(s,s) = C(s,s) - C(s,k)*C(k,s)/C(k,k);
            end
            C = (C+C')/2; 
        else
            break
        end
    end
elseif L == 1
    if s == 0
        v = nu ./ ( prod(M) + nu/tau );
        u = v * h/nu;
        Delta = log(v) + u*conj(u)/v + cnst;
        if Delta>0
            w = u; C = v; s = 1; K = 1;
        end
    else
        Delta = -log(C) - w*conj(w)/C - cnst;
        if Delta>0
            w = 0; C = 0; s = 0; K = 0;
        end
    end
end
end

function [a, theta] = pntFreqEst( eta ,M,D_set_vec, th)
a_theta = @(index,theta)exp(1i*(index-1)*theta');
D = length(M);
d1    = real(eta'*(repmat(a_theta(D_set_vec,th),1,D).*(1i*(D_set_vec-1))))';
d2    = nan(D);
for i = 1:D
    for j = 1:D
        d2(i,j) = real(eta'*(a_theta(D_set_vec,th).*(1i*(D_set_vec(:,i)-1)).*(1i*(D_set_vec(:,j)-1))));
    end
end
e = eig(d2);
dif = prod(e<0);
if dif 
    theta  = th - (d2\d1)';
    kappa =  (Ainv(exp(0.5*diag(inv(d2)))))';
else    
    theta  = th;
    kappa = max(max(abs(eta)))*ones(1,D);
end
a = prod(besseli((D_set_vec-1),repmat(kappa,prod(M),1),1)./besseli(0,repmat(kappa,prod(M),1),1),2).*a_theta(D_set_vec,theta);
end


function [ k ] = Ainv( R )
k   = R; 
in1 = (R<.53); 
in3 = (R>=.85);
in2 = logical(1-in1-in3); 
R1  = R(in1); 
R2  = R(in2); 
R3  = R(in3); 

if ~isempty(R1)
    t      = R1.*R1;
    k(in1) = R1 .* ( 2 + t + 5/6*t.*t );
end
if ~isempty(R2)
    k(in2) = -.4 + 1.39*R2 + 0.43./(1-R2);
end
if ~isempty(R3)
    k(in3) = 1./( R3.*(R3-1).*(R3-3) );
end

end