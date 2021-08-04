%> @file RotationTFYZ.m
%> @brief Implements RotationTFYZ class.
% ==============================================================================
%> @class RotationTFYZ
%> @brief 2-qubit TFYZ rotation gates stored in rotation representation
%>
%> 2-qubit TFYZ rotation gate with circuit representation:
%>
%>      - RX1 -|     |- RX3 -
%>             | RYZ |
%>      - RX2 -|     |- RX4 -
%>
% (C) Copyright Daan Camps and Roel Van Beeumen 2021.  
% ==============================================================================
classdef RotationTFYZ < f3c.qgates.TFTwoAxesRotationGate2
  
  methods
    function [mat] = matrix(obj)
      tmpCirc = qclab.QCircuit( 2 );
      thetas = obj.thetas;
      tmpCirc.push_back( qclab.qgates.RotationX( 0, thetas(1) ) );
      tmpCirc.push_back( qclab.qgates.RotationX( 1, thetas(2) ) );
      tmpCirc.push_back( f3c.qgates.RotationYZ( [0,1], thetas(3), thetas(4) ) );
      tmpCirc.push_back( qclab.qgates.RotationX( 0, thetas(5) ) );
      tmpCirc.push_back( qclab.qgates.RotationX( 1, thetas(6) ) );
      mat = tmpCirc.matrix ;
    end
    
    % toQASM
    function [out] = toQASM(obj, fid, offset)
      if nargin == 2, offset = 0; end
      f3c.IO.qasmRotationTFYZ( fid, obj.qubits + offset, obj.thetas );
      out = 0;
    end
    
    % equalType
    function [bool] = equalType(~, other)
      bool = isa(other,'f3c.qgates.RotationTFYZ');
    end
    
    % label for draw function
    function [label] = label(obj, parameter)
      if nargin < 2, parameter = 'N'; end
      label = 'TFYZ';
      if strcmp(parameter, 'S') % short parameter
        th = obj.thetas;
        label = sprintf([label, '(%.4f,%.4f,%.4f,%.4f,%.4f,%.4f)'], ...
          th(1), th(2), th(3), th(4), th(5), th(6) );
      elseif strcmp(parameter, 'L') % long parameter
        th = obj.thetas;
        label = sprintf([label, '(%.8f,%.8f,%.8f,%.8f,%.8f,%.8f)'], ...
          th(1), th(2), th(3), th(4), th(5), th(6) );
      end
    end
    
  end
  
end % RotationTFYZ