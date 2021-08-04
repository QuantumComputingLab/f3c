%> @file timeStepKitaev.m
%> @brief Generates a Trotter circuit for a single time-step with a Kitaev model.
% ==============================================================================
%> @brief Generates a Trotter circuit for a single time-step with a Kitaev model.
%>
%> @param N       number of spins
%> @param dt      delta t
%> @param Paulis  Pauli string determining the type of gates in the Kitaev
%>                 chain. The length of Paulis should be at least 2 and smaller
%>                 than N. If it is smaller than N-1, the Pauli string is 
%>                 repeated until N-1 is reached. It is a character array 
%>                 containing only 'X',
%>                 'Y', and 'Z'. Valid examples are:
%>                   - 'XY' : XYXYXY...XY Kitaev chain
%>                   - 'YZX' : YZXYZX...YZX Kitaev chain
%>                   - ....
%> @param J  coupling parameters for Kitaev chain. J is either a scalar 
%>            (ordered Hamiltonian) or an array  with N-1 values for the 
%>            coupling strenghts of every nearest neighbor interaction 
%>            (disordered Hamiltonian).
%>
%> @retval timeStepCirc  circuit which implements the Kitaev time-step
%
% (C) Copyright Daan Camps and Roel Van Beeumen 2021.  
% ==============================================================================
function [timeStepCirc] = timeStepKitaev( N, dt, Paulis, J )
  assert( qclab.isNonNegInteger( N ) );
  pLength = length( Paulis );
  assert( pLength > 1 );
  if isscalar( J )
    J = repmat( J, N-1, 1 );
  end
  
  timeStepCirc = qclab.QCircuit( N ) ;
  RXX = @qclab.qgates.RotationXX ;
  RYY = @qclab.qgates.RotationYY ;
  RZZ = @qclab.qgates.RotationZZ ;
  
  % even
  for i = 0:2:N-2
    cPauli = Paulis( mod( i, pLength ) + 1 );
    if strcmp(cPauli, 'X' )
      timeStepCirc.push_back( RXX( [i,i+1], 2*dt*J(i+1) ) );
    elseif strcmp(cPauli, 'Y' )
      timeStepCirc.push_back( RYY( [i,i+1], 2*dt*J(i+1) ) );
    elseif strcmp(cPauli, 'Z' )
      timeStepCirc.push_back( RZZ( [i,i+1], 2*dt*J(i+1) ) );
    else
      warning('unknown Pauli string');
    end
  end
  
  % odd
  for i = 1:2:N-2
    cPauli = Paulis( mod( i, pLength ) + 1 );
    if strcmp(cPauli, 'X' )
      timeStepCirc.push_back( RXX( [i,i+1], 2*dt*J(i+1) ) );
    elseif strcmp(cPauli, 'Y' )
      timeStepCirc.push_back( RYY( [i,i+1], 2*dt*J(i+1) ) );
    elseif strcmp(cPauli, 'Z' )
      timeStepCirc.push_back( RZZ( [i,i+1], 2*dt*J(i+1) ) );
    else
      warning('unknown Pauli string');
    end
  end
end