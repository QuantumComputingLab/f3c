%> @file TFTwoAxesRotationGate2.m
%> @brief Implements TFTwoAxesRotationGate2 class.
% ==============================================================================
%> @class TFTwoAxesRotationGate2
%> @brief Base class for 2-qubit gates for spin Hamiltonians with 
%> nearest-neighbor interactions along two orthogonal axes and a
%> transverse external field:
%>
%>    TFXY, TFXZ, and TFYZ Hamiltonians
%
% (C) Copyright Daan Camps and Roel Van Beeumen 2021.  
% ==============================================================================
classdef TFTwoAxesRotationGate2 < qclab.qgates.QGate2 & qclab.QAdjustable
                        
  properties (Access = protected)
    %> Qubits of this 2-qubit TFTwoAx gate
    qubits_(1,2) int32
    %> left external field rotations
    rotLeft_(1,2) qclab.QRotation
    %> middle coupling rotations
    rotMiddle_(1,2) qclab.QRotation
    %> right external field rotations
    rotRight_(1,2) qclab.QRotation
  end
  
  methods
    % Class constructor  =======================================================
    %> @brief Constructor for 2-qubit rotation gates for 2 axes interactions
    %> with transverse field, i.e.: TFXY, TFXZ, and TFYZ gates. 
    %>
    %> These gates are parametrized by a total of 6 quantum rotations. 
    %> The circuit diagram of a TFTwoAxRotationGate2 looks as follows:
    %>
    %>      --O--|--|--O--
    %>           |  |
    %>      --O--|--|--O--
    %>  
    %> We have two 1-qubit rotations on the left, 2 2-qubit rotations in the
    %> middle, and 2 1-qubit rotations on the right. In TFTwoAxRotationGate2, 
    %> we store these six rotations in, respectively rotLeft_, rotMiddle_, 
    %> and rotRight_.
    %>
    %> The constructor can be used in the following 6 manners:
    %>
    %> 1. TFTwoAxRotationGate2( ) : TFTwoAxRotationGate2 with `qubits` set to
    %>    [0, 1] and all rotations to 0.
    %>
    %> 2. TFTwoAxRotationGate2( qubits ) : TFTwoAxRotationGate2 with `qubits`
    %>    set to `qubits` and all rotations to 0.
    %>
    %> 3. TFTwoAxRotationGate2( qubits, theta1, theta2, theta3, theta4, theta5,
    %>    theta6 ) : TFTwoAxRotationGate2 with `qubits` set to `qubits` and:
    %>
    %>        rotLeft_(1)   = QRotation(theta1), 
    %>        rotLeft_(2)   = QRotation(theta2),
    %>        rotMiddle_(1) = QRotation(theta3), 
    %>        rotMiddle_(2) = QRotation(theta4),
    %>        rotRight_(1)  = QRotation(theta5), 
    %>        rotRight_(2)  = QRotation(theta6),
    %>
    %> 4. TFTwoAxRotationGate2( qubits, theta ) :
    %>    TFTwoAxRotationGate2 with `qubits` set to `qubits` and:
    %>
    %>        rotLeft_(1)   = QRotation(theta(1)), 
    %>        rotLeft_(2)   = QRotation(theta(2)),
    %>        rotMiddle_(1) = QRotation(theta(3)), 
    %>        rotMiddle_(2) = QRotation(theta(4)),
    %>        rotRight_(1)  = QRotation(theta(5)), 
    %>        rotRight_(2)  = QRotation(theta(6)),
    %>
    %> 5. TFTwoAxRotationGate2( qubits, rot1, rot2, rot3, rot4, rot5, rot6 ) :
    %>    TFTwoAxRotationGate2 with `qubits` set to `qubits` and:
    %>
    %>        rotLeft_(1)   = rot1,
    %>        rotLeft_(2)   = rot2,
    %>        rotMiddle_(1) = rot3, 
    %>        rotMiddle_(2) = rot4,
    %>        rotRight_(1)  = rot5, 
    %>        rotRight_(2)  = rot6,
    %>
    %> 6. TFTwoAxRotationGate2( qubits, rot ) :
    %>    TFTwoAxRotationGate2 with `qubits` set to `qubits` and:
    %>
    %>        rotLeft_(1)   = rot(1), 
    %>        rotLeft_(2)   = rot(2),
    %>        rotMiddle_(1) = rot(3), 
    %>        rotMiddle_(2) = rot(4),
    %>        rotRight_(1)  = rot(5), 
    %>        rotRight_(2)  = rot(6),
    %>
    %> 7. TFTwoAxRotationGate2( RotationTFXYMatrix ) :
    %>    TFTwoAxRotationGate2 with `qubits` equal to qubits of RotationTFXYMatrix
    %>    and parameters computed such that the matrix representation agrees
    %>    (in case of TFXY).
    %
    % ==========================================================================
    function obj = TFTwoAxesRotationGate2(qubits, varargin)
      if nargin == 0
        qubits = [0, 1];
      end
      if nargin < 2
        if isa( qubits, 'f3c.qgates.RotationTFXYMatrix' )
          obj.updateByRotationTFXYMatrix( qubits );
          qubits = qubits.qubits;
        else
          obj.rotLeft_(1) = qclab.QRotation( );
          obj.rotLeft_(2) = qclab.QRotation( );
          obj.rotMiddle_(1) = qclab.QRotation( );
          obj.rotMiddle_(2) = qclab.QRotation( );
          obj.rotRight_(1) = qclab.QRotation( );
          obj.rotRight_(2) = qclab.QRotation( );
        end
      elseif isa(varargin{1}, 'qclab.QRotation')
        obj.updateByRotations(varargin{:});
      else
        obj.updateByThetas(varargin{:});
      end
      obj.setQubits( qubits );
    end
    
    %> @brief Update this TFTwoAxesRotationGate2 based on 6 quantum rotations, a
    %> RotationTFXYMatrix or 6 theta values
    function update(obj, varargin)
      assert(~obj.fixed)
      if isa(varargin{1}, 'qclab.QRotation')
        obj.updateByRotations(varargin{:});
      elseif isa(varargin{1}, 'f3c.qgates.RotationTFXYMatrix')
        obj.updateByRotationTFXYMatrix( varargin{1} ) ;
      else
        obj.updateByThetas(varargin{:});
      end
    end
    
    %> @brief Fuse two TFTwoAxesRotationGate2 gates:
    %> side = 'R': -[obj]-[other]-
    %> side = 'L': -[other]-[obj]-
    function fuse(obj, other, side)
      if nargin == 2, side = 'R'; end
      assert(obj.equalType(other) && ~other.fixed && ...
             ~obj.fixed && isequal(obj.qubits,other.qubits));
      objM = f3c.qgates.RotationTFXYMatrix( obj );
      otherM = f3c.qgates.RotationTFXYMatrix( other );
      objM.fuse( otherM, side );
      obj.updateByRotationTFXYMatrix( objM );
    end
    
    %> @brief Matrix multiplication of two TFTwoAxesRotationGate2 gates.
    function [out] = mtimes(obj, other)
      out = copy(obj);
      out.fuse( other );
    end
    
    %> @brief Turnover of three TFTwoAxesRotationGate2 gates
    function [GA, GB, GC] = turnover( G1, G2, G3 )
      assert( G1.equalType(G2) && G1.equalType(G3) );
      G1M = f3c.qgates.RotationTFXYMatrix( G1 );
      G2M = f3c.qgates.RotationTFXYMatrix( G2 );
      G3M = f3c.qgates.RotationTFXYMatrix( G3 );
      [GAM, GBM, GCM] = turnover( G1M, G2M, G3M );
      GA = copy( G2 );
      GA.updateByRotationTFXYMatrix( GAM );
      GB = copy( G1 );
      GB.updateByRotationTFXYMatrix( GBM );
      GC = copy( G2 );
      GC.updateByRotationTFXYMatrix( GCM );
    end
    
    %> @brief rotLeft : left circuit rotations
    function [rotLeft] = rotLeft(obj)
      rotLeft = copy(obj.rotLeft_) ;
    end
    
    function [rotLeftProd] = rotLeftProd(obj) 
      rotLeftProd = obj.rotLeft_(1) * obj.rotLeft_(2) ;
    end
    
    function [rotLeftDiv] = rotLeftDiv(obj) 
      rotLeftDiv = obj.rotLeft_(1) / obj.rotLeft_(2) ;
    end
    
    %> @brief rotMiddle : mtimes and mrdivide of middle circuit rotations
    function [rotMiddle] = rotMiddle(obj)
      rotMiddle = copy(obj.rotMiddle_) ;
    end
    
    %> @brief rotRight : mtimes and mrdivide of right circuit rotations
    function [rotRight] = rotRight(obj)
      rotRight = copy(obj.rotRight_) ;
    end
    
    function [rotRightProd] = rotRightProd(obj)
      rotRightProd = obj.rotRight_(1) * obj.rotRight_(2) ;
    end
    
    function [rotRightDiv] = rotRightDiv(obj) 
      rotRightDiv = obj.rotRight_(1) / obj.rotRight_(2) ;
    end
    
    %> @brief rotations : all 6 circuit rotations ordered from left to right.
    function [rotations] = rotations(obj)
      rotations = [obj.rotLeft, obj.rotMiddle, obj.rotRight];
    end
    
    %> @brief thetasLeft : theta values of left circuit rotations.
    function [thetasLeft] = thetasLeft(obj)
      thetasLeft =  [obj.rotLeft_(1).theta, obj.rotLeft_(2).theta];
    end
    
    %> @brief thetasMiddle : theta values of middle circuit rotations.
    function [thetasMiddle] = thetasMiddle(obj)
      thetasMiddle =  [obj.rotMiddle_(1).theta, obj.rotMiddle_(2).theta];
    end
    
    %> @brief thetasRight : theta values of right circuit rotations.
    function [thetasRight] = thetasRight(obj)
      thetasRight =  [obj.rotRight_(1).theta, obj.rotRight_(2).theta];
    end
    
    %> @brief thetas : all 6 theta values of circuit rotations 
    %> ordered from left to right.
    function [thetas] = thetas(obj)
      thetas = [obj.thetasLeft, obj.thetasMiddle, obj.thetasRight];
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
    
    % equals
    function [bool] = equals(obj,other)
      bool = false;
      if obj.equalType(other)
        bool = isequal( obj.matrix, other.matrix ) ;
      end
    end
    
    %> @brief Checks if `other` is equal to this TFTwoAxRotationGate2.
    function [bool] = eq(obj,other)
      bool = obj.equals(other);
    end
    
    %> @brief Checks if `other` is different from this TFTwoAxRotationGate2.
    function [bool] = ne(obj,other)
      bool = ~eq(obj,other);
    end
    
    % ctranspose
    function objprime = ctranspose( obj )
      objprime = ctranspose@qclab.qgates.QGate2( obj );
      objprime.updateByRotations( inv(obj.rotRight_(1)), inv(obj.rotRight_(2)), ...
                                  inv(obj.rotMiddle_(1)), inv(obj.rotMiddle_(2)), ...
                                  inv(obj.rotLeft_(1)), inv(obj.rotLeft_(2)) );
    end

    %> @brief Convert the TFTwoAxesRotationGate2 to a RotationTFXYMatrix
    function G = f3c.qgates.RotationTFXYMatrix( obj )
      rotMiddleProd = obj.rotMiddle_(1) * obj.rotMiddle_(2) ;
      rotMiddleDiv = obj.rotMiddle_(1) / obj.rotMiddle_(2) ;
      
      theta = obj.rotLeftProd * obj.rotRightProd;
      a = rotMiddleDiv.cos * ( theta.cos - 1i * theta.sin ) ;
      
      theta = obj.rotLeftDiv * obj.rotRightDiv;
      b = rotMiddleProd.cos * (theta.cos - 1i * theta.sin ) ;
      
      theta = obj.rotLeftDiv / obj.rotRightDiv;
      c = -1i*rotMiddleProd.sin * (theta.cos - 1i * theta.sin ) ;
      
      theta = obj.rotLeftProd / obj.rotRightProd;
      d = -1i*rotMiddleDiv.sin * (theta.cos - 1i * theta.sin ) ;
      
      G = f3c.qgates.RotationTFXYMatrix( obj.qubits, [a, b, c, d] );
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
      cp.rotLeft_   = copy(obj.rotLeft_) ;
      cp.rotMiddle_ = copy(obj.rotMiddle_) ;
      cp.rotRight_  = copy(obj.rotRight_) ;
    end
    
    %> @brief Update TFTwoAxRotationGate2 object based on 6 theta values
    %> (in array or not).
    function updateByThetas(obj, varargin)
      assert(nargin == 7 || length(varargin{1}) == 6);
      if nargin == 7
        obj.rotLeft_(1) = qclab.QRotation( varargin{1} );
        obj.rotLeft_(2) = qclab.QRotation( varargin{2} );
        obj.rotMiddle_(1) = qclab.QRotation( varargin{3} );
        obj.rotMiddle_(2) = qclab.QRotation( varargin{4} );
        obj.rotRight_(1) = qclab.QRotation( varargin{5} );
        obj.rotRight_(2) = qclab.QRotation( varargin{6} );
      else
        theta = varargin{1};
        obj.rotLeft_(1)   = qclab.QRotation( theta(1) );
        obj.rotLeft_(2)   = qclab.QRotation( theta(2) );
        obj.rotMiddle_(1) = qclab.QRotation( theta(3) );
        obj.rotMiddle_(2) = qclab.QRotation( theta(4) );
        obj.rotRight_(1)  = qclab.QRotation( theta(5) );
        obj.rotRight_(2)  = qclab.QRotation( theta(6) );
      end
    end
    
    %> @brief Update TFTwoAxRotationGate2 object based on 6 QRotations
    %> (in array or not).
    function updateByRotations(obj,varargin)
      assert(nargin == 7 || length(varargin{1}) == 6);
      if nargin == 7
        obj.rotLeft_(1)   = varargin{1} ;
        obj.rotLeft_(2)   = varargin{2} ;        
        obj.rotMiddle_(1) = varargin{3} ;
        obj.rotMiddle_(2) = varargin{4} ;
        obj.rotRight_(1)  = varargin{5} ;
        obj.rotRight_(2)  = varargin{6} ;
      else
        qrots = varargin{1};
        obj.rotLeft_(1)   = qrots(1) ;
        obj.rotLeft_(2)   = qrots(2) ;        
        obj.rotMiddle_(1) = qrots(3) ;
        obj.rotMiddle_(2) = qrots(4) ;
        obj.rotRight_(1)  = qrots(5) ;
        obj.rotRight_(2)  = qrots(6) ;
      end
    end
    
    %> @brief Update TFTwoAxRotationGate2 object to same parameters as
    %> RotationTFXYMatrix object.
    function updateByRotationTFXYMatrix(obj, RotationTFXYMatrix)
      obj.updateByRotations( RotationTFXYMatrix.rotations );
    end
    
  end
                        
end % TFTwoAxRotationGate2