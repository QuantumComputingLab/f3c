%> @file RotationTFXYMatrix.m
%> @brief Implements RotationTFXYMatrix class.
% ==============================================================================
%> @class RotationTFXYMatrix
%> @brief Class for storing 2-qubit TFXY rotation gates in matrix representation.
%>
%> We use this class to perform fusion and turnover operations for all three
%> types of TFTwoAxes Rotation gates by considering them all as TFXY. We support 
%> transformation between rotation and matrix representations.
%
% (C) Copyright Daan Camps and Roel Van Beeumen 2021.  
% ==============================================================================
classdef RotationTFXYMatrix < qclab.qgates.QGate2 & qclab.QAdjustable
  
  properties (Access = protected)
    %> Qubits of this 2-qubit TFXY gate
    qubits_(1,2) int32
    %> a, b, c, d complex values that represent the TFXY gate
    abcd_(1,4)  double
  end
  
  methods
    % Class constructor  =======================================================
    %> @brief Constructor for 2-qubit TFXY rotation gates in matrix 
    %>  representation.
    %>
    %> We support 3 ways of constructing an object:
    %> 
    %>  1. RotationTFXYMatrix( ) : TFXYMatrix rotation gate with `qubits` set to
    %>     [0,1] and an identity matrix representation.
    %>
    %>  2. RotationTFXYMatrix( qubits ) : TFXYMatrix rotation gate with `qubits` 
    %>     set to `qubits` and an identity matrix representation.
    %>
    %>  3. RotationTFXYMatrix( qubits, abcd ) : TFXYMatrix rotation gate with
    %>     `qubits` set to `qubits` and the matrix element vector to `abcd`.
    %
    % ==========================================================================
    function obj = RotationTFXYMatrix( qubits, abcd )
      if nargin < 1
        qubits = [0, 1];
      end
      if nargin < 2
        abcd = [1, 1, 0, 0];
      end
      obj.setQubits( qubits );
      obj.update( abcd );
    end
    
    %> @brief Update the matrix element vector abcd to `v`
    function update( obj, v )
      obj.abcd_ = v;
    end
    
    % abcd vector
    function [abcd] = abcd(obj)
      abcd = obj.abcd_ ;
    end
    
    % qubit
    function [qubit] = qubit(obj)
      qubit = obj.qubits_(1);
    end
    
    % qubits
    function [qubits] = qubits(obj)
      qubits = obj.qubits_ ;
    end
    
    % setQubits
    function setQubits(obj, qubits)
      assert(qclab.isNonNegIntegerArray(qubits));
      assert(qubits(1) ~= qubits(2));
      obj.qubits_ = sort(qubits(1:2)) ;
    end
    
    function [mat] = matrix(obj)
      mat = [obj.abcd_(1) , 0,            0,        -conj(obj.abcd_(4)) ;
             0,          obj.abcd_(2), -conj(obj.abcd_(3)), 0 ;
             0,          obj.abcd_(3),  conj(obj.abcd_(2)), 0 ;
             obj.abcd_(4),  0,            0,         conj(obj.abcd_(1)) ] ;
    end
    
    % toQASM
    function [out] = toQASM(obj, fid, offset)
      tmp = f3c.qgates.RotationTFXY( obj );
      out = tmp.toQASM( fid, offset );
    end
    
    % equalType
    function [bool] = equalType(~, other)
      bool = isa(other,'f3c.qgates.RotationTFXYMatrix');
    end
    
    % equals
    function [bool] = equals(obj,other)
      bool = false;
      if obj.equalType(other)
        bool = isequal( obj.matrix, other.matrix ) ;
      end
    end
    
    % ctranspose
    function objprime = ctranspose( obj )
      objprime = ctranspose@qclab.qgates.QGate2( obj );
      v = obj.abcd ;
      objprime.update( [conj(v(1)), conj(v(2)), -v(3), -v(4)] );
    end
    
    %> @brief Compute the rotations that make up the circuit for the TFXYMatrix
    %> rotation gate.
    function [rotations] = rotations(obj)
      anglea = angle(obj.abcd_(1)); % -(t1 + t2 + t5 + t6) / 2
      angleb = angle(obj.abcd_(2)); % -(t1 - t2 + t5 - t6) / 2
      anglec = angle(obj.abcd_(3)); % -(t1 - t2 - t5 + t6 + pi) / 2
      angled = angle(obj.abcd_(4)); % -(t1 + t2 - t5 - t6 + pi) / 2
      
      t1p2 = - anglea - angled - pi / 2;
      t1m2 = - angleb - anglec - pi / 2;
      
      t5p6 = - anglea + angled + pi / 2;
      t5m6 = - angleb + anglec + pi / 2;
      
      t3p4 =  atan2( abs(obj.abcd_(3)), abs(obj.abcd_(2)) ) ;
      t3m4 =  atan2( abs(obj.abcd_(4)), abs(obj.abcd_(1)) ) ;
      
      rotations(1) = qclab.QRotation(  ( t1p2 + t1m2 ) / 2 );
      rotations(2) = qclab.QRotation(  ( t1p2 - t1m2 ) / 2 );
      
      rotations(3) = qclab.QRotation( t3p4 + t3m4 );
      rotations(4) = qclab.QRotation( t3p4 - t3m4 );

      rotations(5) = qclab.QRotation( ( t5p6 + t5m6 ) / 2 );
      rotations(6) = qclab.QRotation( ( t5p6 - t5m6 ) / 2 );

    end
    
    %> @brief Fuse two TFXYMatrix rotation gates:
    %> side = 'R': -[obj]-[other]-
    %> side = 'L': -[other]-[obj]-
    function fuse(obj, other, side)
      if nargin == 2, side = 'R'; end
      assert(isa(other,'f3c.qgates.RotationTFXYMatrix') && ~other.fixed && ...
             ~obj.fixed && isequal(obj.qubits,other.qubits));
      if strcmp( side, 'R' )     
        tmp = other.abcd_(1) * obj.abcd_(1) - conj(other.abcd_(4)) * obj.abcd_(4);
        obj.abcd_(4) = other.abcd_(4) * obj.abcd_(1) + ...
                        conj(other.abcd_(1)) * obj.abcd_(4);
        obj.abcd_(1) = tmp;

        tmp = other.abcd_(2) * obj.abcd_(2) - conj(other.abcd_(3)) * obj.abcd_(3);
        obj.abcd_(3) = other.abcd_(3) * obj.abcd_(2) + ...
                        conj(other.abcd_(2)) * obj.abcd_(3);
        obj.abcd_(2) = tmp;
      else % side = 'L'
       tmp = obj.abcd_(1) * other.abcd_(1) - conj(obj.abcd_(4)) * other.abcd_(4);
       obj.abcd_(4) = obj.abcd_(4) * other.abcd_(1) + ...
                        conj(obj.abcd_(1)) * other.abcd_(4);
       obj.abcd_(1) = tmp;
       tmp = obj.abcd_(2) * other.abcd_(2) - conj(obj.abcd_(3)) * other.abcd_(3);
       obj.abcd_(3) = obj.abcd_(3) * other.abcd_(2) + ...
                        conj(obj.abcd_(2)) * other.abcd_(3);
       obj.abcd_(2) = tmp;
      end
    end
    
    %> @brief Matrix multiplication of two RotationTFXYMatrix gates.
    function [out] = mtimes(obj, other)
      out = copy(obj);
      out.fuse( other );
    end
    
    %> @brief Turnover operation on 3 TFXYMatrix rotation gates (implemented in 
    %> different file).
    [GA, GB, GC] = turnover( G1, G2, G3 );
    
    % label for draw function
    function [label] = label(obj, parameter)
      if nargin < 2, parameter = 'N'; end
      label = 'TFXYMat';        
      if strcmp(parameter, 'S') % short parameter
        v = obj.abcd_;
        label = sprintf([label, '(%.4f+%.4fj,%.4f+%.4fj,%.4f+%.4fj,%.4f+%.4fj)'], ...
                real(v(1)), imag(v(1)), real(v(2)), imag(v(2)), ...
                real(v(3)), imag(v(3)), real(v(4)), imag(v(4)));
      elseif strcmp(parameter, 'L') % long parameter
        v = obj.abcd_;
        label = sprintf([label, '(%.8f+%.8fj,%.8f+%.8fj,%.8f+%.8fj,%.8f+%.8fj)'], ...
                real(v(1)), imag(v(1)), real(v(2)), imag(v(2)), ...
                real(v(3)), imag(v(3)), real(v(4)), imag(v(4)));
      end
    end

    %> Convert the RotationTFXYMatrix to a RotationTFXY gate
    function G = f3c.qgates.RotationTFXY( obj )
      G = f3c.qgates.RotationTFXY( obj.qubits, obj.rotations );
    end
    
    %> Convert the RotationTFXYMatrix to a RotationTFXZ gate
    function G = f3c.qgates.RotationTFXZ( obj )
      G = f3c.qgates.RotationTFXZ( obj.qubits, obj.rotations );
    end
  end
  
  methods (Static)
    
    % controlled
    function [bool] = controlled
      bool = false;
    end
    
  end
  
end % RotationTFXYMatrix