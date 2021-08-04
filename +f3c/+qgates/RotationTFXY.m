%> @file RotationTFXY.m
%> @brief Implements RotationTFXY class.
% ==============================================================================
%> @class RotationTFXY
%> @brief 2-qubit TFXY rotation gates stored in rotation representation
%>
%> 2-qubit TFXY rotation gate with circuit representation:
%>
%>      - RZ1 -|     |- RZ3 -
%>             | RXY |
%>      - RZ2 -|     |- RZ4 -
%> 
%> and matrix representation:
%>
%> \f[\begin{bmatrix} a & 0 & 0        & -\bar{d} \\ 
%>                    0 & b & -\bar{c} & 0 \\
%>                    0 & c & \bar{b}  & 0 \\
%>                    b & 0 & 0        & \bar{a} \end{bmatrix},\f]
%> where:
%>
%> \f[a = \cos((\theta_{m,1} - \theta_{m,2})/2) \exp(-i (\theta_{l,1} + 
%>        \theta_{l,2} + \theta_{r,1}) + \theta_{r,2})),\f]
%>
%> \f[b = \cos((\theta_{m,1} + \theta_{m,2})/2) \exp(-i (\theta_{l,1} - 
%>        \theta_{l,2} + \theta_{r,1}) - \theta_{r,2})),\f]
%>
%> \f[c = -i\sin((\theta_{m,1} + \theta_{m,2})/2) \exp(-i (\theta_{l,1} - 
%>        \theta_{l,2} - \theta_{r,1}) + \theta_{r,2})),\f]
%>
%> \f[d = -i\sin((\theta_{m,1} - \theta_{m,2})/2) \exp(-i (\theta_{l,1} + 
%>        \theta_{l,2} - \theta_{r,1}) - \theta_{r,2})),\f]
%>
% (C) Copyright Daan Camps and Roel Van Beeumen 2021.  
% ==============================================================================
classdef RotationTFXY < f3c.qgates.TFTwoAxesRotationGate2
  
  methods

    function [mat] = matrix(obj)           
      tmpCirc = qclab.QCircuit( 2 );
      thetas = obj.thetas;
      tmpCirc.push_back( qclab.qgates.RotationZ( 0, thetas(1) ) );
      tmpCirc.push_back( qclab.qgates.RotationZ( 1, thetas(2) ) );
      tmpCirc.push_back( f3c.qgates.RotationXY( [0,1], thetas(3), thetas(4) ) );
      tmpCirc.push_back( qclab.qgates.RotationZ( 0, thetas(5) ) );
      tmpCirc.push_back( qclab.qgates.RotationZ( 1, thetas(6) ) );
      mat = tmpCirc.matrix ;     
    end
    
    % toQASM
    function [out] = toQASM(obj, fid, offset)
      if nargin == 2, offset = 0; end
      f3c.IO.qasmRotationTFXY( fid, obj.qubits + offset, obj.thetas );
      out = 0;
    end
    
    % equalType
    function [bool] = equalType(~, other)
      bool = isa(other,'f3c.qgates.RotationTFXY');
    end
    
    % label for draw function
    function [label] = label(obj, parameter)
      if nargin < 2, parameter = 'N'; end
      label = 'TFXY';        
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
  
end % RotationTFXY