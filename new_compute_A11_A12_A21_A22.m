function [A11,A12,A21,A22,m] = new_compute_A11_A12_A21_A22(n_par,R,C,Q,tau, r)
    
%% This function computes the A11, A12, A21 and A22 matrices.

%% Stack up the A11 and A12 matrices
A11 = blkdiag(0,-tau(1)); % The A11 matrix.
A12 = [1/Q(1);1/C(1)]; 

for i = 2:n_par
    A12 = blkdiag(A12,[1/Q(i);1/C(i)]);
    A11 = blkdiag(A11,blkdiag(0,-tau(i)));
end

%% % The A22 matrix
% A22 = [R(1)*ones(n_par-1,1),-diag(R(2:n_par)); ones(1,n_par)]; % Computes the A22 matrix

A22 = zeros(n_par, n_par);
% last_row = ones(1, n_par);

% main_diag = -(r(2:n_par) + R(2:n_par));
% main_diag(1) = r(1);

main_diag = -(r + R);
main_diag(1) = r(1);

lower_diag = r(1:n_par-1);
lower_diag(1) = r(1);

A22 = diag(main_diag) + ...
      diag(lower_diag, -1);

% for i = 3:n_par-1
%     A22(i,1) = 0;
% end

for i = 2:n_par
    for j = i+1:n_par
        A22(i, j) = -R(i);
    end
end

A22(1,:) = ones(1,n_par);
disp(A22)

%% Set A21 to be the identity.
A21 = kron(A22(1:n_par-1,:),[0,1]);

%% m is the inv(A22)
m = zeros(n_par,n_par);
Rsum_inv = sum(1./R); 

for ktil = 1:n_par
    for i = 1:n_par-1
        if ktil - i-1 == 0
            m(i+1,i) = (1/R(i+1))^2*(1/Rsum_inv)-1/R(i+1);
        else
            m(ktil,i) = (1/(R(ktil)*R(i+1)))/Rsum_inv;
        end
    end
end

m(n_par,n_par) = 1/(R(n_par)*Rsum_inv);

for i = 1:n_par-1
    m(i,n_par)= 1/(R(i)*Rsum_inv);
end
inv_A22 = inv(A22);
% error = norm(inv_A22-m)

end