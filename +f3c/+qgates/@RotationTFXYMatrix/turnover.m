%> @file turnover.m
%> @brief Turnover operation a pattern of three TFXYMatrix gate objects.
%>
%> Transforms a Vee pattern to a Hat pattern or vice versa. We go from the left
%> pattern to the right pattern or vice versa.
%>
%>@verbatim
%>      G1 G2 G3          G1 G2 G3
%> i   --|-----|--      -----|-----
%> i+1 --|--|--|--  to  --|--|--|--
%> i+2 -----|-----      --|-----|--
%>@endverbatim
%>
%> This function is part of the TFXYMatrix class.
%>
%> @param obj First TFXYMatrix gate in input pattern
%> @param G2 Second TFXYMatrix gate in pattern
%> @param G3 Third TFXYMatrix gate in pattern
%>
%> @retval GA First TFXYMatrix gate in output pattern
%> @retval GB Second TFXYMatrix gate in output pattern
%> @retval GC Third TFXYMatrix gate in output pattern
%
% (C) Copyright Daan Camps and Roel Van Beeumen 2021.  
% ==============================================================================
function [GA, GB, GC] = turnover(obj, G2, G3)
  % Check input types
  assert( isa(G2, 'f3c.qgates.RotationTFXYMatrix') );
  assert( isa(G3, 'f3c.qgates.RotationTFXYMatrix') );
  % Check qubits 
  qubits1 = obj.qubits;
  assert( qubits1(1) + 1 == qubits1(2) );
  qubits2 = G2.qubits;
  assert( qubits2(1) + 1 == qubits2(2) );
  qubits3 = G3.qubits;
  assert( qubits3(1) + 1 == qubits3(2) );
  assert( isequal(qubits1, qubits3) && abs(qubits1(1) - qubits2(1)) == 1 );
  
  % Get abcd vectors
  v1 = obj.abcd_ ;
  v2 = G2.abcd_ ;
  v3 = G3.abcd_ ;

  if qubits2(1) < qubits1(1)
    % Hat to Vee order
    v1(2:3) = conj(v1(2:3));
    v1(3) = -v1(3);
    v2(2:3) = conj(v2(2:3));
    v2(3) = -v2(3);
    v3(2:3) = conj(v3(2:3));
    v3(3) = -v3(3);
  end
  
  % 11 matrix ( & 33 matrix )
  Q11(1,1) =  v3(1) * v2(1)       * v1(1) - conj(v3(4)) * conj(v2(2)) * v1(4) ; 
  Q11(2,1) =  v3(2) * v2(4)       * v1(1) + conj(v3(3)) * conj(v2(3)) * v1(4) ; 
  Q11(1,2) = -v3(1) * conj(v2(4)) * v1(2) - conj(v3(4)) * v2(3)       * v1(3) ;
  Q11(2,2) =  v3(2) * conj(v2(1)) * v1(2) - conj(v3(3)) * v2(2)       * v1(3) ; 
  Q33 = [conj(Q11(2,2)), -conj(Q11(2,1)); -conj(Q11(1,2)), conj(Q11(1,1))];
  
  % 22 matrix ( & 44 matrix )
  Q22(1,1) =  v3(1) * v2(2)       * v1(1) - conj(v3(4)) * conj(v2(1)) * v1(4) ;
  Q22(2,1) =  v3(2) * v2(3)       * v1(1) + conj(v3(3)) * conj(v2(4)) * v1(4) ;
  Q22(1,2) = -v3(1) * conj(v2(3)) * v1(2) - conj(v3(4)) * v2(4)       * v1(3) ;
  Q22(2,2) =  v3(2) * conj(v2(2)) * v1(2) - conj(v3(3)) * v2(1)       * v1(3) ;
  Q44 = [conj(Q22(2,2)), -conj(Q22(2,1)); -conj(Q22(1,2)), conj(Q22(1,1))];
  
  % 41 matrix ( & 23 matrix )
  Q41(1,1) =  v3(3) * v2(4)       * v1(1) - conj(v3(2)) * conj(v2(3)) * v1(4) ;
  Q41(2,1) =  v3(4) * v2(1)       * v1(1) + conj(v3(1)) * conj(v2(2)) * v1(4) ;
  Q41(1,2) =  v3(3) * conj(v2(1)) * v1(2) + conj(v3(2)) * v2(2)       * v1(3) ;
  Q41(2,2) = -v3(4) * conj(v2(4)) * v1(2) + conj(v3(1)) * v2(3)       * v1(3) ;
  Q23 = [conj(Q41(2,2)), -conj(Q41(2,1)); -conj(Q41(1,2)), conj(Q41(1,1))];
  
  % 32 matrix ( & 14 matrix )
  Q32(1,1) =  v3(3) * v2(3)       * v1(1) - conj(v3(2)) * conj(v2(4)) * v1(4) ;
  Q32(2,1) =  v3(4) * v2(2)       * v1(1) + conj(v3(1)) * conj(v2(1)) * v1(4) ;
  Q32(1,2) =  v3(3) * conj(v2(2)) * v1(2) + conj(v3(2)) * v2(1)       * v1(3) ;
  Q32(2,2) = -v3(4) * conj(v2(3)) * v1(2) + conj(v3(1)) * v2(4)       * v1(3) ;
  Q14 = [conj(Q32(2,2)), -conj(Q32(2,1)); -conj(Q32(1,2)), conj(Q32(1,1))];
  
   if norm(Q11) - norm(Q41) < 100 * eps || abs(Q41(2,1)) + abs(Q41(1,2)) == 0 
     
    % diagonalize, anti-diagonalize, diagonalize, anti-diagonalize
    
    % (1) compute U, Y that diagonalize (Q11, Q33)
    [U,Y,D11,~] = f3c.diagonalize22(Q11, Q33);
    
    % (2) compute V to anti-diagonalize (Q41, Q23)
    AD41 = Q41 * Y;
    AD23 = Q23 * Y;
    if abs( AD41(2,1) ) + abs( AD41(2,2) ) > abs( AD23(2,1) ) + abs( AD23(2,2) ) 
      % use AD41
      if abs( AD41(2,2) ) > abs( AD41(2,1) ) %&& AD41(1,2) ~= 0
        V = f3c.rotateToZero( AD41(1,2), AD41(2,2) );
      else
        V = f3c.rotateToZero( AD41(2,1), AD41(1,1) );
        V = rot90( V, 2 );
      end
    else
      % use AD23
      if abs( AD23(2,2) ) > abs( AD23(2,1) )
        V = f3c.rotateToZero( AD23(1,2), AD23(2,2) );
      else
        V = f3c.rotateToZero( AD23(2,1), AD23(1,1) );
        V = rot90( V, 2 );
      end
    end
      
    AD41 = V * AD41;
    
    % (3) compute Z to diagonalize (Q22, Q33)
    D22 = V * Q22;
    D44 = V * Q44;
    % option B
    ZZ = f3c.rotateToZero( D44(2,2), D44(2,1), 'R' );
    DD22 = D22 * ZZ;
    DD44 = D44 * ZZ;
    % option A
    Z = f3c.rotateToZero( D22(2,2), D22(2,1), 'R' );
    D22 = D22 * Z;
    if abs( D11(1,1) - DD22(1,1) ) + abs( D11(2,2) - DD22(2,2) ) < ...
       abs( D11(1,1) - D22(1,1) )  + abs( D11(2,2) - D22(2,2) ) 
      % option B is better than A
      Z = ZZ;
    end
    
  else
    % anti-diagonalize, diagonalize, anti-diagonalize, diagonalize
    
    % (1) compute V, Y that anti-diagonalize ( Q23, Q41 )
    Q23(2,:) = -Q23(2,:);
    Q41(2,:) = -Q41(2,:);
    [V,Y,AD23,AD41] = f3c.diagonalize22( Q23, Q41 );
    V = flipud(V);
    V(:,2) = -V(:,2);
    AD41 = flipud(AD41);
    AD23 = flipud(AD23);
    
    
    % (2) compute U to diagonalize ( Q11, Q33 )
    D11 = Q11 * Y;
    D33 = Q33 * Y;
    if abs( D11(2,1) ) > abs( D33(2,1) )
      U = f3c.rotateToZero( D11(1,1), D11(2,1) );
    else
      U = f3c.rotateToZero( D33(1,1), D33(2,1) );
    end
    D11 = U * D11;
    
    % (3) compute Z to anti-diagonalize ( Q14, Q32 )
    AD32 = U * Q32;
    AD14 = U * Q14;
    % option B
    ZZ = f3c.rotateToZero( AD14(1,2), AD14(1,1), 'R' );
    AAD32 = AD32 * ZZ;
    AAD14 = AD14 * ZZ;
    % option A
    Z = f3c.rotateToZero( AD32(1,2), AD32(1,1), 'R' );
    AD14 = AD14 * Z;
    if abs( AD23(2,1) - AAD14(2,1) ) + abs ( AD23(1,2) - AAD14(1,2) ) < ...
       abs( AD23(2,1) - AD14(2,1) ) + abs ( AD23(1,2) - AD14(1,2) )
       % option B is better than A
       Z = ZZ;
    end
    
  end
  
% Assign outputs
  if qubits2(1) > qubits1(1) % Vee to Hat
    GA = f3c.qgates.RotationTFXYMatrix( qubits2, ...
                     [conj(Y(1,1)), conj(Z(1,1)), conj(Z(1,2)), conj(Y(1,2))] );
    GB = f3c.qgates.RotationTFXYMatrix( qubits1, ...
                                    [D11(1,1), D11(2,2), AD41(1,2), AD41(2,1)]);
    GC = f3c.qgates.RotationTFXYMatrix( qubits2, ...
                     [conj(U(1,1)), conj(V(1,1)), conj(V(1,2)), conj(U(1,2))] );
  else % Hat to Vee
    GA = f3c.qgates.RotationTFXYMatrix( qubits2, ...
                                [conj(Y(1,1)), Z(1,1), -Z(1,2), conj(Y(1,2))] );
    GB = f3c.qgates.RotationTFXYMatrix( qubits1, ...
                      [D11(1,1), conj(D11(2,2)), -conj(AD41(1,2)), AD41(2,1)] );
    GC = f3c.qgates.RotationTFXYMatrix( qubits2, ...
                                [conj(U(1,1)), V(1,1), -V(1,2), conj(U(1,2))] );
  end
end