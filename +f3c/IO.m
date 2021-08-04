%> @file IO.m
%> @brief Implements IO class.
% ==============================================================================
%> @class IO
%> @brief central storage of all QASM strings.
% ==============================================================================
classdef IO
  methods ( Static )
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % QASM I/O
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % qasmRotationXY
    function qasmRotationXY( fid, qubits, thetaX, thetaY )
      fprintf(fid, 'rx(%.15f) q[%d];\n', pi/2, qubits(1) );
      fprintf(fid, 'rx(%.15f) q[%d];\n', pi/2, qubits(2) );
      fprintf(fid, 'cx q[%d], q[%d];\n', qubits(1), qubits(2)); 
      fprintf(fid, 'rx(%.15f) q[%d];\n', thetaX, qubits(1) );
      fprintf(fid, 'rz(%.15f) q[%d];\n', thetaY, qubits(2) );
      fprintf(fid, 'cx q[%d], q[%d];\n', qubits(1), qubits(2)); 
      fprintf(fid, 'rx(%.15f) q[%d];\n', -pi/2, qubits(1) );
      fprintf(fid, 'rx(%.15f) q[%d];\n', -pi/2, qubits(2) );
    end
    % qasmRotationXZ
    function qasmRotationXZ( fid, qubits, thetaX, thetaZ )
      fprintf(fid, 'cx q[%d], q[%d];\n', qubits(1), qubits(2)); 
      fprintf(fid, 'rx(%.15f) q[%d];\n', thetaX, qubits(1) );
      fprintf(fid, 'rz(%.15f) q[%d];\n', thetaZ, qubits(2) );
      fprintf(fid, 'cx q[%d], q[%d];\n', qubits(1), qubits(2)); 
    end
    
    % qasmRotationYZ
    function qasmRotationYZ( fid, qubits, thetaY, thetaZ )
      fprintf(fid, 'rz(%.15f) q[%d];\n', pi/2, qubits(1) );
      fprintf(fid, 'rz(%.15f) q[%d];\n', pi/2, qubits(2) );
      fprintf(fid, 'cx q[%d], q[%d];\n', qubits(1), qubits(2)); 
      fprintf(fid, 'rx(%.15f) q[%d];\n', thetaY, qubits(1) );
      fprintf(fid, 'rz(%.15f) q[%d];\n', thetaZ, qubits(2) );
      fprintf(fid, 'cx q[%d], q[%d];\n', qubits(1), qubits(2)); 
      fprintf(fid, 'rz(%.15f) q[%d];\n', -pi/2, qubits(1) );
      fprintf(fid, 'rz(%.15f) q[%d];\n', -pi/2, qubits(2) );
    end
    
    % qasmRotationTFXY
    function qasmRotationTFXY( fid, qubits, thetas )
      fprintf(fid, 'rz(%.15f) q[%d];\n', thetas(1), qubits(1) );
      fprintf(fid, 'rz(%.15f) q[%d];\n', thetas(2), qubits(2) );
      fprintf(fid, 'rx(%.15f) q[%d];\n', pi/2, qubits(1) );
      fprintf(fid, 'rx(%.15f) q[%d];\n', pi/2, qubits(2) );
      fprintf(fid, 'cx q[%d], q[%d];\n', qubits(1), qubits(2)); 
      fprintf(fid, 'rx(%.15f) q[%d];\n', thetas(3), qubits(1) );
      fprintf(fid, 'rz(%.15f) q[%d];\n', thetas(4), qubits(2) );
      fprintf(fid, 'cx q[%d], q[%d];\n', qubits(1), qubits(2)); 
      fprintf(fid, 'rx(%.15f) q[%d];\n', -pi/2, qubits(1) );
      fprintf(fid, 'rx(%.15f) q[%d];\n', -pi/2, qubits(2) );
      fprintf(fid, 'rz(%.15f) q[%d];\n', thetas(5), qubits(1) );
      fprintf(fid, 'rz(%.15f) q[%d];\n', thetas(6), qubits(2) );
    end
    
    % qasmRotationTFXZ
    function qasmRotationTFXZ( fid, qubits, thetas )
      fprintf(fid, 'ry(%.15f) q[%d];\n', thetas(1), qubits(1) );
      fprintf(fid, 'ry(%.15f) q[%d];\n', thetas(2), qubits(2) );
      fprintf(fid, 'cx q[%d], q[%d];\n', qubits(1), qubits(2)); 
      fprintf(fid, 'rx(%.15f) q[%d];\n', thetas(3), qubits(1) );
      fprintf(fid, 'rz(%.15f) q[%d];\n', thetas(4), qubits(2) );
      fprintf(fid, 'cx q[%d], q[%d];\n', qubits(1), qubits(2)); 
      fprintf(fid, 'ry(%.15f) q[%d];\n', thetas(5), qubits(1) );
      fprintf(fid, 'ry(%.15f) q[%d];\n', thetas(6), qubits(2) );
    end
    
    % qasmRotationTFYZ
    function qasmRotationTFYZ( fid, qubits, thetas )
      fprintf(fid, 'rx(%.15f) q[%d];\n', thetas(1), qubits(1) );
      fprintf(fid, 'rx(%.15f) q[%d];\n', thetas(2), qubits(2) );
      fprintf(fid, 'rz(%.15f) q[%d];\n', pi/2, qubits(1) );
      fprintf(fid, 'rz(%.15f) q[%d];\n', pi/2, qubits(2) );
      fprintf(fid, 'cx q[%d], q[%d];\n', qubits(1), qubits(2)); 
      fprintf(fid, 'rx(%.15f) q[%d];\n', thetas(3), qubits(1) );
      fprintf(fid, 'rz(%.15f) q[%d];\n', thetas(4), qubits(2) );
      fprintf(fid, 'cx q[%d], q[%d];\n', qubits(1), qubits(2)); 
      fprintf(fid, 'rz(%.15f) q[%d];\n', -pi/2, qubits(1) );
      fprintf(fid, 'rz(%.15f) q[%d];\n', -pi/2, qubits(2) );
      fprintf(fid, 'rx(%.15f) q[%d];\n', thetas(5), qubits(1) );
      fprintf(fid, 'rx(%.15f) q[%d];\n', thetas(6), qubits(2) );
    end
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  end
end