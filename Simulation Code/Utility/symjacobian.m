function G = symjacobian(y,x)
% Compute Symbolic Jacobian that can be evaluated as an LFR
    G = jacobian(y,x);
    G = char(G);
    G = strrep(G,'], [',';'); 
    G = strrep(G,'matrix([[','['); 
    G = strrep(G,']])',']');
end