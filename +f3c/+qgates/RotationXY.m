%> @file RotationXY.m
%> @brief Implements RotationXY class.
% ==============================================================================
%> @class RotationXY
%> @brief 2-qubit XY rotation gates stored in rotation representation
%>
%> 2-qubit XY rotation gate with matrix representation:
%>
%> \f[\begin{bmatrix} c_{-}    & 0        & 0        & -i s_{-} \\ 
%>                    0        & c_{+}    & -i s_{+} & 0 \\
%>                    0        & -i s_{+} & c_{+}    & 0 \\
%>                    -i s_{-} & 0        & 0        & c_{-} \end{bmatrix},\f]
%> where:
%>
%> \f[c_{-} = \cos((\theta_{1} - \theta_{2})/2),\f]
%>
%> \f[c_{+} = \cos((\theta_{1} + \theta_{2})/2),\f]
%>
%> \f[s_{-} = \sin((\theta_{1} - \theta_{2})/2),\f]
%>
%> \f[s_{+} = \sin((\theta_{1} + \theta_{2})/2),\f]
%>
% (C) Copyright Daan Camps and Roel Van Beeumen 2021.  
% ==============================================================================
classdef RotationXY < f3c.qgates.TwoAxesRotationGate2
  
  methods
    % matrix
    function [mat] = matrix(obj)
      tplus = obj.rotations_(1) * obj.rotations_(2) ;
      tminus = obj.rotations_(1) / obj.rotations_(2) ;
      
      mat = [ tminus.cos,     0,             0,             -1i*tminus.sin ;
              0,              tplus.cos,     -1i*tplus.sin, 0 ;
              0,              -1i*tplus.sin, tplus.cos,     0 ;
              -1i*tminus.sin, 0,             0,             tminus.cos ];
    end
    
    % toQASM
    function [out] = toQASM(obj,fid, offset)
      if nargin == 2, offset = 0; end
      f3c.IO.qasmRotationXY( fid, obj.qubits + offset, ...
                               obj.thetas(1), obj.thetas(2) );
      out = 0;
    end
    
    % equals
    function [bool] = equals(obj,other)
      bool = false;
      if obj.equalType(other)
        bool = isequal( obj.matrix, other.matrix ) ;
      end
    end
    
    % equalType
    function [bool] = equalType(~, other)
      bool = isa(other,'f3c.qgates.RotationXY');
    end
    
    % label for draw function
    function [label] = label(obj, parameter)
      if nargin < 2, parameter = 'N'; end
      label = 'XY';        
      if strcmp(parameter, 'S') % short parameter
        th = obj.thetas;
        label = sprintf([label, '(%.4f,%.4f)'], th(1), th(2) );
      elseif strcmp(parameter, 'L') % long parameter
        th = obj.thetas;
        label = sprintf([label, '(%.8f,%.8f)'], th(1), th(2) );
      end
    end
    
  end
end % RotationXY

