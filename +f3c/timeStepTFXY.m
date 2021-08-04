%> @file timeStepTFXY.m
%> @brief Generates a Trotter circuit for a single time-step with a TFXY model.
% ==============================================================================
%> @brief Generates a Trotter circuit for a single time-step with a TFXY model.
%>
%> @param N   number of spins
%> @param dt  delta t
%> @param Jx  coupling along x-axis. Jx is either a scalar (ordered) or an array
%>            with N-1 values for the coupling strenghts of every nearest
%>            neighbor interaction (disordered).
%> @param Jy  coupling along y-axis. Jy is either a scalar (ordered) or an array
%>            with N-1 values for the coupling strenghts of every nearest
%>            neighbor interaction (disordered).
%> @param hz  transverse field strength along z-axis. hz is either a scalar
%>            (ordered) or an array with N-1 values for the external field
%>            strengths at every spin (disordered).
%>
%> @retval timeStepMat  circuit which implements the time-step where the TFXY
%>                      gates are stored as matrices.
%> @retval timeStepRot  circuit which implements the time-step where the TFXY
%>                      gates are stored as angles.
%
% (C) Copyright Daan Camps and Roel Van Beeumen 2021.  
% ==============================================================================
function [timeStepMat, timeStepRot] = timeStepTFXY( N, dt, Jx, Jy, hz )
  assert( qclab.isNonNegInteger( N ) );
  if isscalar( Jx )
    Jx = repmat( Jx, N-1, 1 );
  end
  assert( length(Jx) == N-1 );
  if isscalar( Jy )
    Jy = repmat( Jy, N-1, 1 );
  end
  assert( length(Jy) == N-1 );
  if isscalar( hz )
    hz = repmat( hz, N, 1 );
  end
  assert( length(hz) == N );
  
  timeStepMat = qclab.QCircuit( N );
  timeStepRot = qclab.QCircuit( N );
  TFXY = @f3c.qgates.RotationTFXY;
  TFXYM = @f3c.qgates.RotationTFXYMatrix;
  
  % even
  for i = 0:2:N-2
    gate = TFXY( [i,i+1], 2*dt*hz(i+1), 2*dt*hz(i+2), 2*dt*Jx(i+1), ...
                                                      2*dt*Jy(i+1), 0, 0 );
    timeStepRot.push_back( gate );
    gate = TFXYM( gate );
    timeStepMat.push_back( gate );
  end
  
  % odd
  for i = 1:2:N-3
    gate = TFXY( [i,i+1], 0, 0, 2*dt*Jx(i+1), 2*dt*Jy(i+1), 0, 0 );
    timeStepRot.push_back( gate );
    gate = TFXYM( gate );
    timeStepMat.push_back( gate );
  end
  
  if mod(N,2) == 1 % odd number of gates
    gate = TFXY( [N-2,N-1], 0, 2*dt*hz(N), 2*dt*Jx(N-1), 2*dt*Jy(N-1), 0, 0 );
    timeStepRot.push_back( gate );
    gate = TFXYM( gate );
    timeStepMat.push_back( gate );
  end
  
end