%> @file diagonalize22.m
%> @brief Compute a diagonalization for a 2 x 2 diagonalizable matrix pencil.
%>
%> The function is typically used for diagonalizable pencils of the form:
%> %> \f[\begin{bmatrix} a   &  b\\ c & d \end{bmatrix} - \lambda 
%>    \begin{bmatrix} \bar{d} & -\bar{c} \\ -\bar{b} & \bar{a} \end{bmatrix},\f]
%>
%> @param A 2 x 2 complex matrix
%> @param B 2 x 2 complex matrix
%>
%> @retval Q 2 x 2 complex matrix containing left eigenvectors
%> @retval Z 2 x 2 complex matrix containing right eigenvectors
%> @retval DA 2 x 2 complex diagonal matrix: A in eigenbasis
%> @retval DB 2 x 2 complex diagonal matrix: B in eigenbasis
%
% (C) Copyright Daan Camps and Roel Van Beeumen 2021.  
% ==============================================================================
function [Q,Z,A,B] = diagonalize22(A,B)
  maxit = 100;
  stepit = 3;
  
  tolA = 12 * eps * norm( A );
  tolB = 12 * eps * norm( B );
  
  Q = eye(2); Z = eye(2);
  for iter=1:maxit
    % compute eigenvalues to get the shift.
    if mod(iter, stepit) == 0
      % exceptional shift to try to get out of stagnation
      if A(2,2) ~= 0
        a1 = A(2,2); b = B(2,2); a2 = 0;
      else
        a1 = A(1,1); b = B(1,1); a2 = 0;
      end
    else
      [a1, a2, b] = f3c.eig22(A, B);
    end
    
    if abs(a1) > abs(a2)
      H = b * A - a1 * B;
    else
      H = b * A - a2 * B;
    end
  
    if abs( H(1,1) ) + abs( H(2,1) ) >  abs( H(1,2) ) + abs( H(2,2) )
      QQ = f3c.rotateToZero( H(1,1), H(2,1) );
    else
      QQ = f3c.rotateToZero( H(1,2), H(2,2) );
    end
  
    A = QQ * A;
    B = QQ * B;
  
    if abs( A(2,1) ) + abs( A(2,2) ) >  abs( B(2,1) ) + abs( B(2,2) ) 
      ZZ = f3c.rotateToZero( A(2,2), A(2,1), 'R' );
    else
      ZZ = f3c.rotateToZero( B(2,2), B(2,1), 'R' );
    end
  
    A = A * ZZ;
    B = B * ZZ;
    
    Q = QQ * Q;
    Z = Z * ZZ;
    
    % check exit condition
    if (abs(A(2,1)) + abs(A(1,2)) < tolA && ...
        abs(B(2,1)) + abs(B(1,2)) < tolB )
        break
    end
    
    % method did not converge
    if iter == maxit
      warning(['maximum number of iterations in diagonalize22 reached, ', ...
               'results may be inaccurate.\n', ...
               '  |A(1,2)| + |A(2,1)| = %.4e  (tolA = %.4e)\n', ...
               '  |B(1,2)| + |B(2,1)| = %.4e  (tolB = %.4e)\n'], ...
               abs(A(2,1)) + abs(A(1,2)), tolA,  abs(B(2,1)) + abs(B(1,2)), tolB);
    end
  end
end


