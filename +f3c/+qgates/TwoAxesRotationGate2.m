%> @file TwoAxesRotationGate2.m
%> @brief Implements TwoAxesRotationGate2 class.
% ==============================================================================
%> @class TwoAxesRotationGate2
%> @brief Base class for 2-qubit gates for spin Hamiltonians with 
%> nearest-neighbor interactions along two orthogonal axes:
%>
%>    XY, XZ, and YZ Hamiltonians
%
% (C) Copyright Daan Camps and Roel Van Beeumen 2021.  
% ==============================================================================
classdef TwoAxesRotationGate2 < qclab.qgates.QGate2 & qclab.QAdjustable
  
  properties (Access = protected)
    %> Qubits of this 2-qubit TwoAxRotation gate
    qubits_(1,2) int32
    %> two rotations of this TwoAx gate
    rotations_(1,2) qclab.QRotation
  end
  
  methods
    % Class constructor  =======================================================
    %> @brief Constructor for 2-qubit rotation gates for 2 axes interactions,
    %> i.e.: XY, XZ, and YZ gates. 
    %>
    %> These gates are parametrized by two quantum rotations. 
    %> The circuit diagram of a TwoAxGate2 looks as follows:
    %>
    %>      --|--|--
    %>        |  |
    %>      --|--|--
    %>
    %> The constructor can be used in the following 6 manners:
    %>
    %> 1. TwoAxesRotationGate2( ) : TwoAxesRotationGate2 with `qubits` set to
    %>    [0, 1] and all rotations to 0.
    %>
    %> 2. TwoAxesRotationGate2( qubits ) : TwoAxesRotationGate2 with `qubits` 
    %>    set to `qubits` and  the rotations to 0.
    %>
    %> 3. TwoAxesRotationGate2( qubits, theta1, theta2 ) : TwoAxesRotationGate2
    %>    with `qubits` set to `qubits` and rotations to \f$\theta_1/2\f$ and
    %>    \f$\theta_2/2\f.
    %>
    %> 4. TwoAxRotationGate2( qubits, rot1, rot2 ) : TwoAxesRotationGate2 with
    %>    `qubits` set to `qubits` and rotations to `rot1` and `rot2`.
    %
    % ==========================================================================
    function obj = TwoAxesRotationGate2(qubits, varargin)
      if nargin == 0
        qubits = [0, 1];
      end
      obj.setQubits( qubits );
      if nargin == 3
        if isa( varargin{1}, 'qclab.QRotation')
          obj.rotations_(1) = varargin{1} ;
          obj.rotations_(2) = varargin{2} ;
        else
          obj.rotations_(1) = qclab.QRotation( varargin{1} );
          obj.rotations_(2) = qclab.QRotation( varargin{2} );
        end
      else
        obj.rotations_(1) = qclab.QRotation( );
        obj.rotations_(2) = qclab.QRotation( );
      end
    end
    
    %> @brief Update this TwoAxesRotationGate2
    function update(obj, var1, var2)
      assert(~obj.fixed);
      if isa(var1, 'qclab.QRotation')
        obj.rotations_(1) = var1 ;
        obj.rotations_(2) = var2 ;
      else
        obj.rotations_(1) = qclab.QRotation( var1 );
        obj.rotations_(2) = qclab.QRotation( var2 );
      end
    end
    
    %> @brief Fuse two TwoAxesRotationGate2 gates. Other is positioned on `side`
    %> of obj.
    function fuse(obj, other, side)
      if nargin == 2, side = 'R'; end
      assert(obj.equalType(other) && ~other.fixed && ...
             ~obj.fixed && isequal(obj.qubits,other.qubits));
      % fuse two rotations
      obj.rotations_(1) = obj.rotations_(1) * other.rotations_(1) ;
      obj.rotations_(2) = obj.rotations_(2) * other.rotations_(2) ;
    end
    
    %> @brief Turnover of pattern of three same-type TwoAxesRotationGate2 gates.
    function [GA, GB, GC] = turnover(obj, G2, G3)
      % Check qubits and type
      assert( obj.equalType(G2) && obj.equalType(G3) );
      qubits1 = obj.qubits;
      assert( qubits1(1) + 1 == qubits1(2) );
      qubits2 = G2.qubits;
      assert( qubits2(1) + 1 == qubits2(2) );
      qubits3 = G3.qubits;
      assert( qubits3(1) + 1 == qubits3(2) );
      assert( isequal(qubits1, qubits3) && abs(qubits1(1) - qubits2(1)) == 1 );

      % Create output objects of correct type
      GA = feval( class(obj) );
      GB = feval( class(obj) );
      GC = feval( class(obj) );
      
      % two SU(2) turnovers
      [RA2, RB1, RC2] = ...
       qclab.turnoverSU2( obj.rotations_(1), G2.rotations_(2), G3.rotations_(1) );
      [ RA1, RB2, RC1 ] = ...
       qclab.turnoverSU2( obj.rotations_(2), G2.rotations_(1), G3.rotations_(2) );
      
      % update rotations
      GA.update( RA1, RA2 );
      GB.update( RB1, RB2 );
      GC.update( RC1, RC2 );
      
      % update qubits
      GA.setQubits( qubits2 );
      GB.setQubits( qubits1 );
      GC.setQubits( qubits2 );
    end
    
    %> @brief Matrix multiplication of two TwoAxesRotationGate2 gates.
    function [out] = mtimes(obj, other)
      out = copy(obj);
      out.fuse( other );
    end
    
    %> @brief Two theta values of TwoAxesRotationGate2
    function [thetas] = thetas(obj)
      thetas = [ obj.rotations_(1).theta, obj.rotations_(2).theta ];
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
    
    %> @brief Checks if `other` is equal to this TwoAxesRotationGate2.
    function [bool] = eq(obj,other)
      bool = obj.equals(other);
    end
    
    %> @brief Checks if `other` is different from this TwoAxesRotationGate2.
    function [bool] = ne(obj,other)
      bool = ~eq(obj,other);
    end
    
    % ctranspose
    function objprime = ctranspose( obj )
      objprime = ctranspose@qclab.qgates.QGate2( obj );
      objprime.update( inv(obj.rotations_(1)), inv(obj.rotations_(2)) );
    end
    
  end
  
   methods (Static)
    
    % controlled
    function [bool] = controlled
      bool = false;
    end
    
  end
  
  methods (Abstract)
    
    %> checks if other is of the same type as this
    [bool] = equalType(obj, other)
    
  end
  
  methods ( Access = protected )
    
    function cp = copyElement(obj)
      cp = copyElement@matlab.mixin.Copyable( obj );
      cp.rotations_ = copy(obj.rotations_) ;
    end
    
  end
  
end % TwoAxesRotationGate2