function [ HPsi ] = propagate_1D_HSB(input,T,V,Psi)    
    
    % apply the full system-bath hamiltonian
    
    % T, V, and number operator
    HPsi = ifft(bsxfun(@times,T,fft(Psi))) + bsxfun(@times,V,Psi) + bsxfun(@times,input.epsilon,Psi);
    
    % apply the mixing part of the bath
    for a = 1:2^input.bathN
        HPsi(:,a) = HPsi(:,a) + sum(bsxfun(@times,input.Je_coupling,Psi(:,input.relax_index(a,:))),2);
    end
    
    % apply the dephasing part of the path (the full dephasing matrix Okl
    % is 0 at the boundaries)
    for a = 2:2^input.bathN-1
        ind_dephase = find(input.dephase_index(a,:));
        HPsi(:,a) = HPsi(:,a) + sum(bsxfun(@times,input.cij(a,ind_dephase),Psi(:,input.dephase_index(a,ind_dephase))),2);
    end
end