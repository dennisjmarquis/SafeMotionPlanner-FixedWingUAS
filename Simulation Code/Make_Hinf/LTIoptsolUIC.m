function [R, S, scal, gamma] = LTIoptsolUIC(A, B1, B2, C1, C2, D11, D12, D21, Lambda, N_UIC, gamma)

nx = length(A(:,:,1));
nz = size(C1,1);
nw = size(D21,2);

Ns = null([C2 D21]);
U1 = Ns(1:nx,:);
U2 = Ns(nx+1:nx+nw,:);

Nr = null([B2' D12']);
V1 = Nr(1:nx,:);
V2 = Nr(nx+1:nx+nz,:);

% defing the decision variables (symmetric real matrix)
for kk = 1:N_UIC+1
Rvar{kk} = sdpvar(nx);
Svar{kk} = sdpvar(nx);
end
evar = sdpvar(1);
f1var = sdpvar(1);
f2var = sdpvar(1);
pvar = sdpvar(1);
tvar = sdpvar(1);

if isempty(gamma) % mimimization
    gvar = sdpvar(1);
else   % no minimization
    gvar = gamma;
end

eps = 1e-9;
L = [A B1;C1 D11];

%% Define LMIs

LMIsys = [];

LMIsys = LMIsys + ...
        [evar+f1var+f2var+eps <= 2*gvar];

LMIsys = LMIsys + ...
        [Lambda'*Svar{1}*Lambda + eps*eye(size(Lambda,2)) <= f1var*eye(size(Lambda,2))];
% LMIsys = LMIsys + ...
%         [Svar{1} + eps*eye(size(Svar,1)) <= f1var*eye(size(Svar,1))];

LMIsys = LMIsys + ...
        [[pvar 1; 1 f2var] >= 0];

LMIsys = LMIsys + ...
        [[tvar 1; 1 evar] >= 0];

for kk = 1:N_UIC
    LMIsys = LMIsys + ...
        [[V1;V2]'*(L*[Rvar{kk} zeros(nx,nw);zeros(nw,nx) pvar*eye(nw)]*L' - [Rvar{kk+1} zeros(nx,nz);zeros(nz,nx) evar*eye(nz)])*[V1;V2] + eps*eye(size(Nr,2)) <= 0];

    LMIsys = LMIsys + ...
        [[U1;U2]'*(L'*[Svar{kk+1} zeros(nx,nz);zeros(nz,nx) tvar*eye(nz)]*L - [Svar{kk} zeros(nx,nw);zeros(nw,nx) f2var*eye(nw)])*[U1;U2] + eps*eye(size(Ns,2)) <= 0];

    LMIsys = LMIsys + ...
        [[Rvar{kk} eye(nx);eye(nx) Svar{kk}]>=0];
end

       LMIsys = LMIsys + ...
        [[V1;V2]'*(L*[Rvar{N_UIC+1} zeros(nx,nw);zeros(nw,nx) pvar*eye(nw)]*L' - [Rvar{N_UIC+1} zeros(nx,nz);zeros(nz,nx) evar*eye(nz)])*[V1;V2] + eps*eye(size(Nr,2)) <= 0];

    LMIsys = LMIsys + ...
        [[U1;U2]'*(L'*[Svar{N_UIC+1} zeros(nx,nz);zeros(nz,nx) tvar*eye(nz)]*L - [Svar{N_UIC+1} zeros(nx,nw);zeros(nw,nx) f2var*eye(nw)])*[U1;U2] + eps*eye(size(Ns,2)) <= 0];

    LMIsys = LMIsys + ...
        [[Rvar{N_UIC+1} eye(nx);eye(nx) Svar{N_UIC+1}]>=0];


% ops = sdpsettings('showprogress',1,'solver','sdpt3','sdpt3.maxit',maxit,'shift',1e-8);
ops = sdpsettings('verbose',1,'showprogress',0,'solver','mosek');
solvesdp(LMIsys,gvar,ops);

for kk=1:N_UIC+1
    R{kk} = double(Rvar{kk});
    S{kk} = double(Svar{kk});
end
gamma = double(gvar);
scal.e = double(evar);
scal.f1 = double(f1var);
scal.f2 = double(f2var);
scal.p = double(pvar);
scal.t = double(tvar);

end