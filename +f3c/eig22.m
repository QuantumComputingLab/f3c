%> @file eig22.m
%> @brief Compute the eigenvalues of a 2 x 2 matrix pencil.
%>
%> The eigenvalues are  \f$\lambda_1 = \alpha_1/\beta, 
%> \lambda_2 = \alpha_2/\beta\f$
%> 
%> @param A  2 x 2 complex matrix
%> @param B  2 x 2 complex matrix
%> 
%> @retval alpha1 numerator of first generalized eigenvalue
%> @retval alpha2 numerator of second generalized eigenvalue
%> @retval beta denominator of both generalized eigenvalues
%
% (C) Copyright Daan Camps and Roel Van Beeumen 2021.  
% ==============================================================================
function [alpha1, alpha2, beta] = eig22(A, B)
  
  AB11 = A(1,1) * B(2,2) - A(2,1) * B(1,2);
  AB12 = A(1,2) * B(2,2) - A(2,2) * B(1,2);
  AB21 = A(2,1) * B(1,1) - A(1,1) * B(2,1);
  AB22 = A(2,2) * B(1,1) - A(1,2) * B(2,1);
  
  t = (AB11 + AB22) / 2;
  d = sqrt(t^2 + AB12 * AB21 - AB11 * AB22);
  
  beta = B(1,1) * B(2,2) - B(1,2) * B(2,1);
  
  alpha1 = t + d;
  alpha2 = t - d;
  
end