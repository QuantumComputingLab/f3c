%> @file TriangleCircuit.m
%> @brief Implements TriangleCircuit class.
% ==============================================================================
%> @class TriangleCircuit
%> @brief Class for representing a triangle quantum circuit for Hamiltonian
%>  simulation
%>
%> In a TriangleCircuit, the gates appear in the following pattern:
%>
%>           |
%>          | |
%>         | | |
%>        | | | |
%>       | | | | |
%>
%> We support two possible orderings of the gates.
%>
%>    1) Ascending ordering:
%> 
%> 
%>              05                  /
%>            04  09               / /
%>          03  08  12            / / /
%>        02  07  11  14         / / / /
%>      01  06  10  13  15      / / / / /  
%>
%>
%>    2) Descending ordering:
%>
%>
%>              11                \
%>            07  12             \ \
%>          04  08  13          \ \ \
%>        02  05  09  14       \ \ \ \
%>      01  03  06  10  15    \ \ \ \ \
%>
% (C) Copyright Daan Camps and Roel Van Beeumen 2021.  
% ==============================================================================
classdef TriangleCircuit < qclab.QCircuit
  
  properties (Access = protected)
    %> flag that indicates how the TriangleCircuit is stored: ascending or
    %descening order.
    ascend_(1,1)    logical
  end
  
  methods
    
    % Constructor ==============================================================
    %> Constructor for SquareCircuit objects
    %>
    %>   1. TriangleCircuit() : create a 1-qubit, empty TriangleCircuit
    %>
    %>   2. TriangleCircuit( circuit, ascend ) : create a TriangleCircuit 
    %>      from a regular circuit of correct dimensions. If ascend is `true`,
    %>      circuit should be a triangle with ascending ordering; otherwise with
    %>      descending ordering.
    %
    % ==========================================================================
    function obj = TriangleCircuit( circuit, ascend ) 
      obj = obj@qclab.QCircuit( 1, 0 );
      if nargin == 0
        obj.ascend_ = true ;
      elseif nargin == 2
        obj.update( circuit, ascend );
      end
    end
    
    %> @brief Update this TriangleCircuit with a QCircuit of right dimensions.
    function update( obj, circuit, ascend )
      assert( isa( circuit, 'qclab.QCircuit' ) );
      nbQubits = double( circuit.nbQubits );
      nbGates = circuit.nbGates ;
      assert( (nbQubits * (nbQubits - 1))/2 == nbGates );
      % assert that all the gates are ordered as expected and of compatible
      % type.
      types = cell( nbQubits - 1, 1 );
      if ascend % circuit should be ascend triangle
        % l_0 contains gates on all qubits, store for verification
        for q = 0:nbQubits-2
          [linearIndex] = obj.computeAscendLinearIndex(nbQubits, q, 0);
          gateq = circuit.gateHandles( linearIndex );
          assert( isa( gateq, 'qclab.qgates.QGate2' ) );
          qubitsq = gateq.qubits ;
          assert( qubitsq(1) == q );
          assert( qubitsq(1) + 1 == qubitsq(2) );
          types{q+1} = class(gateq);
        end
        % check all other layers
        for l = 1:nbQubits - 2
          for q = l:nbQubits - 2
            [linearIndex] = obj.computeAscendLinearIndex(nbQubits, q, l);
            gateql = circuit.gateHandles( linearIndex );
            qubitsql = gateql.qubits ;
            assert( qubitsql(1) == q );
            assert( qubitsql(1) + 1 == qubitsql(2) );
            assert( isa( gateql, types{q + 1} ) );
          end
        end
        % all OK for ascending
      else % circuit should be descending triangle
        % l_nbQubits-2 contains gates on all qubits, store for verification
        for q = 0:nbQubits-2
          [linearIndex] = obj.computeDescendLinearIndex(nbQubits, q, nbQubits-2);
          gateq = circuit.gateHandles( linearIndex );
          assert( isa( gateq, 'qclab.qgates.QGate2' ) );
          qubitsq = gateq.qubits ;
          assert( qubitsq(1) == q );
          assert( qubitsq(1) + 1 == qubitsq(2) );
          types{q+1} = class(gateq);
        end
        % check all other layers
        for l = 0:nbQubits - 3
          for q = nbQubits - l - 2:nbQubits - 2
            [linearIndex] = obj.computeDescendLinearIndex(nbQubits, q, l);
            gateql = circuit.gateHandles( linearIndex );
            qubitsql = gateql.qubits ;
            assert( qubitsql(1) == q );
            assert( qubitsql(1) + 1 == qubitsql(2) );
            assert( isa( gateql, types{q + 1} ) );
          end
        end
        % all OK for descending
      end      
      % assign properties
      obj.nbQubits_ = nbQubits ;
      obj.ascend_ = ascend ;
      obj.gates_ = circuit.gateHandles ;
    end
    
    %> @brief Returns if this TriangleCircuit is stored in ascending order.
    function [bool] = ascend( obj )
      bool = obj.ascend_ ;
    end
    
    %> @brief Returns if this TriangleCircuit is stored in descending order.
    function [bool] = descend( obj )
      bool = ~obj.ascend_ ;
    end
    
    %> @brief Reorders a descending TriangleCircuit to an ascending
    %> TriangleCircuit. Involves a copy of all gates. No effect if already
    %> ascending.
    function makeAscend( obj )
      if ~obj.ascend_
        gatesCopy = copy( obj.gates_ );
        j = 1;
        for l = 0:obj.nbQubits - 2
          for q = obj.nbQubits - 2:-1:l
            [linearIndex] = obj.ascendLinearIndex( q, l );
            obj.gates_( j ) = gatesCopy( linearIndex );
            j = j + 1;
          end
        end
        obj.ascend_ = true;
      end
    end
    
    
    %> @brief Reorders an ascending TriangleCircuit to descending
    %> TriangleCircuit. Involves a copy of all gates. No effect if already
    %> descending.
    function makeDescend( obj )
      if obj.ascend_
        gatesCopy = copy( obj.gates_ );
        j = 1;
        for l = 0:obj.nbQubits - 2
          for q = obj.nbQubits - l - 2:obj.nbQubits - 2
            [linearIndex] = obj.descendLinearIndex( q, l );
            obj.gates_( j ) = gatesCopy( linearIndex );
            j = j + 1;
          end
        end
        obj.ascend_ = false;
      end
    end
    
    
    %> @brief Compute the linear index of a gate from a qubit and layer index in 
    %> a TriangleCircuit where the qubit and layer index are based on an 
    %> ascending ordering of the triangle.
    %>
    %> @param obj   TriangleCircuit
    %> @param qubitIndex  qubit index of gate ranging from layerIndex to
    %>                    obj.nbQubits-2
    %> @param layerIndex layer index of gate ranging from 0 to obj.nbQubits-2
    %>
    %> @retval linearIndex linear index of the gate in the TriangleCircuit.
    function [linearIndex] = ascendLinearIndex(obj, qubitIndex, layerIndex)
      nbQubits = double( obj.nbQubits );
      if obj.ascend_ 
        linearIndex = obj.computeAscendLinearIndex( nbQubits, qubitIndex, ...
                                                      layerIndex );
      else
        % modify input ascending layerIndex to descending layerIndex
        layerIndex = layerIndex - qubitIndex + nbQubits - 2;
        linearIndex = obj.computeDescendLinearIndex( nbQubits, qubitIndex, ...
                                                       layerIndex );
      end
      
    end
    
    %> @brief Compute the linear index of a gate from a qubit and layer index 
    %> in a TriangleCircuit where the qubit and layer index are based on a
    %> descending ordering of the triangle.
    %>
    %> @param obj   TriangleCircuit
    %> @param qubitIndex  qubit index of gate ranging from obj.nbQubits -
    %>                    layerIndex - 2 to obj.nbQubits - 2
    %> @param layerIndex  layer index of gate ranging from 0 to obj.nbQubits-2
    %>
    %> @retval linearIndex linear index of the specified gate in the
    %>                      TriangleCircuit.
    function [linearIndex] = descendLinearIndex(obj, qubitIndex, layerIndex)
      nbQubits = double( obj.nbQubits );
      if obj.ascend_
        % modify input descending layerIndex to ascending layerIndex
        layerIndex = layerIndex + qubitIndex - nbQubits + 2;
        linearIndex = obj.computeAscendLinearIndex( nbQubits, qubitIndex, ...
                                                      layerIndex );
      else
        linearIndex = obj.computeDescendLinearIndex( nbQubits, qubitIndex, ...
                                                       layerIndex );
      end
    end
    
    %> @brief Compute the qubit and layer index of a gate in a TriangleCircuit
    %> from the linear index. The qubit and layer indices are compatible with
    %> the ordering of the TriangleCircuit.
    %>
    %> @param obj TriangleCircuit
    %> @param linearIndex linear gate index for TriangleCircuit between
    %>                    1:obj.nbGates
    %> 
    %> @retval qubitIndex qubit index of gate ranging from layerIndex to
    %>                    obj.nbQubits-2
    %> @retval layerIndex layer index of gate ranging from 0 to obj.nbQubits-2
    %>    according to the ordering of the TriangleCircuit (ascend or descend ).
    function [qubitIndex, layerIndex] = QLIndex(obj, linearIndex)
      nbQubits = double( obj.nbQubits );
      if obj.ascend_
        [qubitIndex, layerIndex] = obj.computeAscendQLIndex(nbQubits,linearIndex);
      else
        [qubitIndex, layerIndex] = obj.computeDescendQLIndex(nbQubits, linearIndex);
      end  
    end
    
    %> @brief merge gate `gate` from side `side` with this TriangleCircuit.
    %>
    %> A copy of the gate is create such that it is not modified in this
    %> routine, only the triangle circuit is modified.
    %>
    %> @param obj TriangleCircuit
    %> @param gate gate that can be merged with TriangleCircuit
    %> @param side 'L' or 'R' if the gate is positioned on the left or right
    %>        side of the triangle. The default is 'R'.
    function merge( obj, gate, side )
      gate = copy( gate ); 
      qubit = gate.qubit;
      assert( qubit < obj.nbQubits - 1 );
      linearIndex = obj.ascendLinearIndex( qubit, 0 );
      assert( isa( gate, class(obj.gates_( linearIndex ) ) ) );
      if nargin < 3, side = 'R'; end
      if strcmp( side, 'L' ) 
        currentLayer = 0 ; % gate starts at layer zero
        for q = qubit:obj.nbQubits - 3
          linearIndex2 = obj.ascendLinearIndex( q + 1, currentLayer ); 
          linearIndex3 = obj.ascendLinearIndex( q, currentLayer );
          [ obj.gates_(linearIndex2), obj.gates_(linearIndex3), gate ] = ...
            turnover( gate, obj.gates_(linearIndex2), obj.gates_(linearIndex3) );
          currentLayer = currentLayer + 1;
        end
        % fuse from left
        linearIndex = obj.ascendLinearIndex( obj.nbQubits - 2, currentLayer ); 
        obj.gates_(linearIndex).fuse( gate, 'L' );
        
      else % Right side
        currentLayer = obj.nbQubits - 2; % gate starts at rightmost layer
        for q = qubit:obj.nbQubits - 3
          LinearIndex1 = obj.descendLinearIndex( q, currentLayer );
          linearIndex2 = obj.descendLinearIndex( q + 1, currentLayer );
          [ gate, obj.gates_(LinearIndex1), obj.gates_(linearIndex2) ] = ...
            turnover( obj.gates_(LinearIndex1), obj.gates_(linearIndex2), gate );
          currentLayer = currentLayer - 1;
        end
        % fuse from right
        linearIndex = obj.descendLinearIndex( obj.nbQubits - 2, currentLayer ); 
        obj.gates_(linearIndex).fuse( gate, 'R' );
      end
    end
    
    %> @brief Replace a gate at `pos` with a gate of equal type.
    function replace( obj, pos, gate )
      assert( qclab.isNonNegInteger( pos - 1 ) );
      assert( pos <= obj.nbGates );
      assert( length(gate) == 1 );
      assert( gate.equalType( obj.gates_( pos ) ) );
      obj.gates_( pos ) = gate ;
    end
    
    %> @brief Conversion to a qclab QCircuit `C` containing a copy of `obj`.
    function C = qclab.QCircuit( obj )
      C = qclab.QCircuit( obj.nbQubits );
      C.insert( 1:nbGates, obj.gates );
    end
    
    %> @brief Conversion to a SquareCircuit `SqC`. Doesn't modify `obj`.
    function SqC = f3c.SquareCircuit( obj )
      % Triangle to Square algorithm based is based on ascending Triangle.
      
      nbQubits = double( obj.nbQubits );
      C = qclab.QCircuit( nbQubits );
      gatesNew = repmat( qclab.qgates.Identity, 1, obj.nbGates );
      gatesCopy = copy( obj.gates_ );
      
      if mod( nbQubits, 2 ) == 0 % even qubits
        turnoverLayerRng = nbQubits-2:-2:2;
        moveLayerRng = 1:2:nbQubits-3;
      else % odd qubits
        turnoverLayerRng = nbQubits-2:-2:1;
        moveLayerRng = 0:2:nbQubits-3;
      end
      
      % Turnover the odd slanted layers from right to left through Triangle
      for turnoverLayer = turnoverLayerRng
        for qubit = nbQubits - 2:-1:turnoverLayer
         % get the gate from the Triangle
         [linearIndex] = obj.ascendLinearIndex( qubit, turnoverLayer );
         gate = gatesCopy( linearIndex );
         initQubit = gate.qubit ;
         for layer = turnoverLayer - 1:-1:0
           linearIndex1 = obj.ascendLinearIndex( gate.qubit, layer );
           linearIndex2 = obj.ascendLinearIndex( gate.qubit - 1, layer );
           [gate, gatesCopy(linearIndex1), gatesCopy(linearIndex2)] = ...
             turnover(gatesCopy(linearIndex1), gatesCopy(linearIndex2), gate);
         end
         SqIndex = f3c.SquareCircuit.computeLinearIndex( nbQubits, ...
                                      gate.qubit, nbQubits - initQubit - 2);
         gatesNew( SqIndex ) = gate;
        end
      end
      
      % Add the even slanted layers without turnovers
      if mod( nbQubits, 2 ) == 0 % even qubits
        % Move slanted layer0
        for qubit = nbQubits - 2:-1:0
          TrIndex = obj.ascendLinearIndex( qubit, 0 );
          SqIndex = f3c.SquareCircuit.computeLinearIndex( nbQubits, ...
                         qubit, nbQubits - qubit - 2 );
          gatesNew( SqIndex ) = gatesCopy( TrIndex );
        end
      end
      
      for moveLayer = moveLayerRng
        for qubit = nbQubits - 2:-1:moveLayer
          TrIndex = obj.ascendLinearIndex( qubit, moveLayer );
          SqIndex = f3c.SquareCircuit.computeLinearIndex( nbQubits, ...
                          qubit, nbQubits - qubit + moveLayer - 1 );
          gatesNew( SqIndex ) = gatesCopy( TrIndex );
        end
      end
      
      % assign output
      C.insert( 1:length(gatesNew), gatesNew );
      SqC = f3c.SquareCircuit( C );
    end   
    
  end
  
  methods (Static)
    % Compute circuit indices
    
    %> @brief Compute the linear index of a gate from a qubit and layer index 
    %> for an ascending TriangleCircuit.
    %>
    %> @param nbQubits   number of qubits of ascending TriangleCircuit
    %> @param qubitIndex  qubit index of gate ranging from layerIndex to
    %>                    nbQubits-2
    %> @param layerIndex layer index of gate ranging from 0 to nbQubits-2
    %>
    %> @retval linearIndex linear index of the gate in the TriangleCircuit.
    function [linearIndex] = computeAscendLinearIndex(nbQubits, qubitIndex, ...
                                                       layerIndex)
      assert( layerIndex < nbQubits - 1 );
      assert( qubitIndex >= layerIndex && qubitIndex < nbQubits - 1 );
      linearIndex = ( layerIndex * (2*nbQubits - layerIndex - 1) ) / 2 + ...
                       nbQubits - qubitIndex - 1;        
    end
    
    %> @brief Compute the linear index of a gate from a qubit and layer index 
    %> for a descending TriangleCircuit
    %>
    %> @param nbQubits   number of qubits of descending TriangleCircuit
    %> @param qubitIndex  qubit index of gate ranging from nbQubits -
    %>                    layerIndex - 2 to nbQubits - 2
    %> @param layerIndex  layer index of gate ranging from 0 to nbQubits-2
    %>
    %> @retval linearIndex linear index of the specified gate in the
    %>                      TriangleCircuit.
    function [linearIndex] = computeDescendLinearIndex(nbQubits, qubitIndex, ...
                                                        layerIndex)
      assert( layerIndex < nbQubits - 1 );
      assert( qubitIndex + layerIndex >= nbQubits - 2 && ...
               qubitIndex < nbQubits - 1 );
      linearIndex = (layerIndex * (layerIndex + 1 ) ) / 2 + ...
                        qubitIndex - nbQubits + layerIndex + 3;
    end
    
    %> @brief Compute the qubit and layer index of a gate from a linear index
    %> for an ascending TriangleCircuit
    %>
    %> @param nbQubits   number of qubits of ascending TriangleCircuit
    %> @param linearIndex linear gate index for TriangleCircuit between
    %>                    1:(nbQubits * (nbQubits - 1))/2
    %> 
    %> @retval qubitIndex qubit index of gate ranging from layerIndex to
    %>                    nbQubits-2
    %> @retval layerIndex layer index of gate ranging from 0 to nbQubits-2
    function [qubitIndex, layerIndex] = computeAscendQLIndex(nbQubits, ...
                                                              linearIndex)
      assert( linearIndex <= nbQubits * (nbQubits - 1) / 2);
      % layerIndex
      layerIndex = nbQubits - ...
        (sqrt( (2 * nbQubits - 1)^2 - 8 * linearIndex ) + 1 ) / 2;
      r = layerIndex - floor(layerIndex);
      layerIndex = floor(layerIndex);
      if r == 0, layerIndex = layerIndex - 1; end
      % qubitIndex
      qubitIndex =  ( layerIndex * (2*nbQubits - layerIndex - 1) ) / 2 + ...
        nbQubits - linearIndex - 1;
    end
    
    %> @brief Compute the qubit and layer index of a gate from a linear index
    %> for a descending TriangleCircuit
    %>
    %> @param nbQubits   number of qubits of descending TriangleCircuit
    %> @param linearIndex linear gate index for TriangleCircuit between
    %>                    1:(nbQubits * (nbQubits - 1))/2
    %>
    %> @retval qubitIndex qubit index of gate ranging from nbQubits -
    %>                    layerIndex - 2 to nbQubits - 2
    %> @retval layerIndex layer index of gate ranging from 0 to nbQubits-2
    function [qubitIndex, layerIndex] = computeDescendQLIndex(nbQubits, ...
                                                               linearIndex)
      assert( linearIndex <= nbQubits * (nbQubits - 1) / 2);
      % layerIndex
      layerIndex = (sqrt( 1 + 8 * linearIndex ) - 1 ) / 2;
      r = layerIndex - floor(layerIndex);
      layerIndex = floor(layerIndex);
      if r == 0, layerIndex = layerIndex - 1; end
      % qubitIndex
      qubitIndex = linearIndex - (layerIndex * (layerIndex + 1 ) ) / 2 + ...
                    nbQubits - layerIndex - 3;
    end
    
    
    % Disable QCircuit functionality to limit modifiability of TriangleCircuits
    
    function setOffset(obj, offset)
      assert( false );
    end
    
    function [gates] = gateHandles( obj, pos )
      assert( false );
    end
    
    %
    % Override QCircuit Modifiers for TriangleCircuit
    %
    
    function clear(obj)
      assert( false );
    end
    
    function erase(obj)
      assert( false );
    end
    
    function insert(obj, pos, gates)
      assert( false );
    end
    
    function push_back( obj, gate )
      assert( false );
    end
    
    function pop_back( obj )
      assert( false );
    end
    
    
  end
end % TriangleCircuit

