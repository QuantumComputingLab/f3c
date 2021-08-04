%> @file IsingCircuit.m
%> @brief Implements IsingCircuit class.
% ==============================================================================
%> @class IsingCircuit
%> @brief Useful for compressing circuits for classical Ising models
%> with one and two qubit gates of the same type to depth O(1)
% ==============================================================================
classdef IsingCircuit < qclab.QObject & qclab.QAdjustable
  properties (Access = protected )
    nbQubits_(1,1)                  int32
    type_(1,1)                      char
    circuitRotationGates1_          qclab.QCircuit
    circuitRotationGates2_          qclab.QCircuit
  end
  
  methods
    % 
    function obj = IsingCircuit( nbQubits, type )
      assert(qclab.isNonNegInteger(nbQubits-1)) ;
      assert( type == 'X' || type == 'Y' || type == 'Z' );
      obj.nbQubits_ = nbQubits ;
      obj.type_ = type ;
      % initialize circuits
      obj.circuitRotationGates1_ = qclab.QCircuit( nbQubits, 0 );
      obj.circuitRotationGates2_ = qclab.QCircuit( nbQubits, 0 );
      % add gates with zero angles
      if strcmp( obj.type_, 'X' )
        OneQRot = @qclab.qgates.RotationX ;
        TwoQRot = @qclab.qgates.RotationXX ;
      elseif strcmp( obj.type_, 'Y' )
        OneQRot = @qclab.qgates.RotationY ;
        TwoQRot = @qclab.qgates.RotationYY ;
      else
        OneQRot = @qclab.qgates.RotationZ ;
        TwoQRot = @qclab.qgates.RotationZZ ;
      end
      for i = 0:nbQubits-1
        obj.circuitRotationGates1_.push_back( OneQRot( i ) );
      end
      for i = 0:2:nbQubits-2
          obj.circuitRotationGates2_.push_back( TwoQRot( [i, i+1] ) );
      end
      for i = 1:2:nbQubits-2
          obj.circuitRotationGates2_.push_back( TwoQRot( [i, i+1] ) );
      end
    end
    
    % type: X, Y, Z
    function [type] = type(obj)
      type = obj.type_;
    end
    
    % change type: X, Y, Z
    function setType(obj, type )
      assert( type == 'X' || type == 'Y' || type == 'Z' );
      if type ~= obj.type_
         if strcmp( type, 'X' )
           OneQRot = @qclab.qgates.RotationX ;
           TwoQRot = @qclab.qgates.RotationXX ;
         elseif strcmp( type, 'Y' )
           OneQRot = @qclab.qgates.RotationY ;
           TwoQRot = @qclab.qgates.RotationYY ;
         else
           OneQRot = @qclab.qgates.RotationZ ;
           TwoQRot = @qclab.qgates.RotationZZ ;
         end
         for i = 1:obj.nbQubits
           gate = obj.circuitRotationGates1_.gates( i );
           obj.circuitRotationGates1_.replace( i, OneQRot( gate ) );
         end
         for i = 1:obj.nbQubits-1
           gate = obj.circuitRotationGates2_.gates( i );
           obj.circuitRotationGates2_.replace( i, TwoQRot( gate ) );
         end
         obj.type_ = type;
      end
    end
    
    % merge one qubit rotation gates
    function  mergeOneQRot( obj, gate )
      if strcmp( obj.type_, 'X' )
        assert( isa( gate, 'qclab.qgates.RotationX') );
      elseif strcmp( obj.type_, 'Y' )
        assert( isa( gate, 'qclab.qgates.RotationY') );
      else
        assert( isa( gate, 'qclab.qgates.RotationZ') );
      end
      qubit = gate.qubit ;
      assert( qubit < obj.nbQubits_ );
      obj.circuitRotationGates1_.gateHandles( qubit + 1 ).fuse( gate );
    end
    
    % merge two qubit rotation gates
    function mergeTwoQRot( obj, gate )
      if strcmp( obj.type_, 'X' )
        assert( isa( gate, 'qclab.qgates.RotationXX') );
      elseif strcmp( obj.type_, 'Y' )
        assert( isa( gate, 'qclab.qgates.RotationYY') );
      else
        assert( isa( gate, 'qclab.qgates.RotationZZ') );
      end
      qubit = gate.qubit ;
      assert( qubit < obj.nbQubits_ - 1 );
      if mod( qubit, 2 ) == 0
        index = ( qubit / 2 ) + 1;
      else
        index = ( obj.nbQubits_+ qubit ) / 2;
      end
      obj.circuitRotationGates2_.gateHandles( index ).fuse( gate );
    end
    
    % nbQubits
    function [nbQubits] = nbQubits(obj)
      nbQubits = obj.nbQubits_;
    end
    
    % qubits
    function [qubits] = qubits(obj)
      qubits = 0:obj.nbQubits_-1;
    end
    
    % matrix
    function [mat] = matrix(obj)
      mat = obj.circuitRotationGates1_.matrix * ...
              obj.circuitRotationGates2_.matrix ;
    end
    
    % apply
    function [mat] = apply(obj, side, op, nbQubits, mat, offset )
      assert( nbQubits >= obj.nbQubits_ );
      if nargin == 5, offset = 0; end
      mat = obj.circuitRotationGates1_.apply( side, op, nbQubits, mat, offset );
      mat = obj.circuitRotationGates2_.apply( side, op, nbQubits, mat, offset );
    end
    
    % toQASM
    function [out] = toQASM(obj, fid, offset)
      if nargin == 2, offset = 0; end
      out = obj.circuitRotationGates1_.toQASM( fid, offset );
      if ( out ~= 0 ), return; end
      out = obj.circuitRotationGates2_.toQASM( fid, offset );
    end
    
    % equals
    function [bool] = equals(obj, other)
      bool = false;
      if isa(other, 'f3c.IsingCircuit')
        bool = isequal( other.matrix, obj.matrix );
      end
    end
    
    % draw
    function [out] = draw(obj, fid, parameter, offset)
      if nargin < 2, fid = 1; end
      if nargin < 3, parameter = 'N'; end
      if nargin < 4, offset = 0; end
      qclab.drawCommands ; % load draw commands
      circuitCell1 = obj.circuitRotationGates1_.draw( 0, parameter, 0 );
      circuitCell2 = obj.circuitRotationGates2_.draw( 0, parameter, 0 );
      circuitCell = cell(3*obj.nbQubits,1);
      
      for i = 1:3*obj.nbQubits
        circuitCell{ i } = strcat( circuitCell1{ i }, circuitCell2{ i } );
      end
      
      if fid > 0
        qclab.drawCellArray( fid, circuitCell, obj.qubits + offset );
        out = 0;
      else
        out = circuitCell ;
      end
      
    end
  end
  
  methods (Static)
    
    % controlled
    function [bool] = controlled
      bool = false;
    end
    
    % qubit
    function [qubit] = qubit(obj)
      qubit = 0;
    end
    
    % setQubit
    function setQubit(~)
      assert( false );
    end
    
    % setQubits
    function setQubits(~)
      assert( false );
    end
  end
  
  methods ( Access = protected )
    
    %> @brief Override copyElement function to allow for correct deep copy of
    %> handles.
    function cp = copyElement(obj)
      cp = copyElement@matlab.mixin.Copyable( obj );
      cp.circuitRotationGates1_ = copy( obj.circuitRotationGates1_ );
      cp.circuitRotationGates2_ = copy( obj.circuitRotationGates2_ );
    end
    
  end
  
end % class IsingCircuit