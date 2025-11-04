function [Ak, Bk, Ck, Dk] = LTIcontrollerUIC(R, S, scal, A, B1, B2, C1, C2, D11, D12, D21, gamma)

    nx = length(A(:,:,1));
    nw = size(D11,2);
    nu = size(D12,2);
    nz = size(C1,1);
    ny = size(C2,1);
    
    N = length(R);
    
    R{N+1} = R{N};
    S{N+1} = S{N};
    
    e = scal.e;
    f1 = scal.f1;
    f2 = scal.f2;
    p = scal.p;
    t = scal.t;


% ... Set up matrices for LMI we nee to solve H+Q'J'P+P'JQ < 0 ...
   P = [zeros(nx) eye(nx) zeros(nx,2*nx) zeros(nx,nw) zeros(nx,nz);...
       B2' zeros(nu,nx) zeros(nu,2*nx) zeros(nu,nw) e^(-1/2)*D12'];

   Q = [zeros(nx,2*nx) zeros(nx) eye(nx) zeros(nx,nw) zeros(nx,nz);...
       zeros(ny,2*nx) C2 zeros(ny,nx) f2^(-1/2)*D21 zeros(ny,nz)];

   for kk = 1:N+1
       E{kk} = (R{kk}-inv(S{kk}))^(1/2);
   end
    
   for kk = 1:N
       H{kk} = [-R{kk+1} -E{kk+1} A zeros(nx) f2^(-1/2)*B1 zeros(nx,nz);...
           -E{kk+1}' -eye(nx) zeros(nx) zeros(nx) zeros(nx,nw) zeros(nx,nz);...
           A' zeros(nx) -S{kk} S{kk}*E{kk} zeros(nx,nw) e^(-1/2)*C1';...
           zeros(nx) zeros(nx) E{kk}'*S{kk} -eye(nx)-E{kk}'*S{kk}*E{kk} zeros(nx,nw) zeros(nx,nz);...
           f2^(-1/2)*B1' zeros(nw,nx) zeros(nw,nx) zeros(nw,nx) -eye(nw) (e*f2)^(-1/2)*D11';...
           zeros(nz,nx) zeros(nz,nx) e^(-1/2)*C1 zeros(nz,nx) (e*f2)^(-1/2)*D11 -eye(nz)];
       Jvar{kk} = sdpvar(size(P,1),size(Q,1));
    end
    
    
    fprintf('Solving Controller LMI...\n')
    
    eps=1e-9;
    
    LMIsys = [];
    
    for kk=1:N
    LMIsys = LMIsys + [H{kk}+Q'*Jvar{kk}'*P+P'*Jvar{kk}*Q + eps*eye(size(H{kk},1)) <= 0];
    end
    
    ops = sdpsettings('verbose',0,'solver','mosek');
    diagnostics = optimize(LMIsys,norm(Jvar{N},2),ops);
    
    for kk=1:N
        J=double(Jvar{kk});
        Ak{kk} = J(1:nx,1:nx);
        Bk{kk} = J(1:nx,nx+1:nx+ny);
        Ck{kk} = J(nx+1:nx+nu,1:nx);
        Dk{kk} = J(nx+1:nx+nu,nx+1:nx+ny);
    end

end
