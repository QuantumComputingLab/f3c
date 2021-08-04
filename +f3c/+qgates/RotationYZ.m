%> @file RotationYZ.m
%> @brief Implements RotationYZ class.
% ==============================================================================
%> @class RotationYZ
%> @brief 2-qubit YZ rotation gates stored in rotation representation
%>
%> 2-qubit YZ rotation gate with matrix representation:
%>
%> \f[\begin{bmatrix} c_{-}    & 0        & 0        & s_{-} \\ 
%>                    0        & c_{+}    & s_{+} & 0 \\
%>                    0        & s_{+} & c_{+}    & 0 \\
%>                    s_{-} & 0        & 0        & c_{-} \end{bmatrix},\f]
%> where:
%>
%> \f[c_{-} = e^{-i \theta_2 / 2} \cos(\theta_{1}/2),\f]
%>
%> \f[c_{+} = e^{i \theta_2 / 2} \cos(\theta_{1}/2),\f]
%>
%> \f[s_{-} = i e^{-i \theta_2 / 2} \sin(\theta_{1}/2),\f]
%>
%> \f[s_{+} = -i e^{i \theta_2 / 2} \sin(\theta_{1}/2),\f]
%>
% (C) Copyright Daan Camps and Roel Van Beeumen 2021.  
% ==============================================================================
classdef RotationYZ < f3c.qgates.TwoAxesRotationGate2
  
  methods
    % matrix
    function [mat] = matrix(obj)
      cminus = (obj.rotations_(2).cos - 1i * obj.rotations_(2).sin) * ...
                 obj.rotations_(1).cos;
      cplus =  (obj.rotations_(2).cos + 1i * obj.rotations_(2).sin) * ...
                 obj.rotations_(1).cos;
      sminus = (1i * obj.rotations_(2).cos + obj.rotations_(2).sin) * ...
                 obj.rotations_(1).sin;         
      splus = (-1i * obj.rotations_(2).cos + obj.rotations_(2).sin) * ...
                 obj.rotations_(1).sin;  
      mat = [ cminus, 0,     0,     sminus ;
              0,      cplus, splus, 0 ;
              0,      splus, cplus, 0 ;
              sminus, 0,     0,     cminus ];
    end
    
    % toQASM
    function [out] = toQASM(obj,fid, offset)
      if nargin == 2, offset = 0; end
      f3c.IO.qasmRotationYZ( fid, obj.qubits + offset, ...
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
      bool = isa(other,'f3c.qgates.RotationYZ');
    end
    
    % label for draw function
    function [label] = label(obj, parameter)
      if nargin < 2, parameter = 'N'; end
      label = 'YZ';        
      if strcmp(parameter, 'S') % short parameter
        th = obj.thetas;
        label = sprintf([label, '(%.4f,%.4f)'], th(1), th(2) );
      elseif strcmp(parameter, 'L') % long parameter
        th = obj.thetas;
        label = sprintf([label, '(%.8f,%.8f)'], th(1), th(2) );
      end
    end
    
  end
end % RotationYZ

