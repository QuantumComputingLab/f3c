%> @file TFIMTriangleCircuit.m
%> @brief Implements TFIMTriangleCircuit class.
% ==============================================================================
%> @class TFIMTriangleCircuit
%>
%> In a TFIMTriangleCircuit, the one- and two-spin gates appear in the following
%> pattern:
%>
%> Ascending ordering:
%>
%>            | O         |
%>        | O |       | O |         | O       |
%>    | O |       | O |         | O |     | O |     | O     |
%>  O |         O |           O |       O |       O |     O |   O
%>
%>
%>
%> Descending ordering:
%>
%>
%>                                        |             O |
%>                    |       O |         | O |           | O |
%>      |     O |     | O |     | O |         | O |           | O |
%>  O   | O     | O       | O       | O           | O             | O
%>
%>
% (C) Copyright Daan Camps and Roel Van Beeumen 2021.  
% ==============================================================================
classdef TFIMTriangleCircuit < qclab.QCircuit
  properties (Access = protected)
    %> flag that indicates how the TFIMTriangleCircuit is stored: ascending or
    %descening order.
    ascend_(1,1)    logical
  end
  
  methods
    function obj = TFIMTriangleCircuit( circuit, ascend )
      obj = obj@qclab.QCircuit(1, 0 );
      if nargin == 0
        obj.ascend_ = true;
        obj.push_back( qclab.qgates.RotationZ(0, 0 ) );
      else
        obj.update( circuit, ascend );
      end
    end
    
    function update( obj, circuit, ascend )
      assert( isa( circuit, 'qclab.QCircuit' ) );
      nbQubits = double( circuit.nbQubits );
      nbGates = circuit.nbGates;
      assert( nbQubits * (2*nbQubits - 1) == nbGates );
      type1 = class( circuit.gateHandles(1) );
      type2 = class( circuit.gateHandles(2) );
      if ascend
        i = 1;
        for l = 0:2*nbQubits-2
          q = nbQubits - 1;
          for j = 0:2*nbQubits-2-l
            gatei = circuit.gateHandles(i);
            assert( gatei.qubit == q );
            i = i + 1;
            if mod(j,2) == 0
              assert( isa( gatei, type1 ) );
              q = q - 1;
            else
              assert( isa( gatei, type2 ) );
            end
          end
        end % -- all OK for ascending
        obj.ascend_ = true;
      else
        i = 1;
        for l = 0:2*nbQubits-2
          if mod(l,2) == 0
            q = nbQubits - 1 - (l/2);
            for j = 0:l
              gatei = circuit.gateHandles(i);
              assert( gatei.qubit == q );
              i = i + 1;
              if mod(j,2) == 0
                assert( isa( gatei, type1 ) );
              else
                assert( isa( gatei, type2 ) );
                q = q + 1;
              end
            end
          else
            q = nbQubits - 1 - ((l+1)/2);
            for j=0:l
              gatei = circuit.gateHandles(i);
              assert( gatei.qubit == q );
              i = i + 1;
              if mod(j,2) == 0
                assert( isa( gatei, type2 ) );
                q = q + 1;
              else
                assert( isa( gatei, type1 ) );
              end
            end
          end
        end % -- all OK for descending
        obj.ascend_ = false;
      end
      obj.nbQubits_ = nbQubits ;
      obj.gates_ = circuit.gateHandles() ;
    end
    
    %> @brief Returns if this TFIMTriangleCircuit is stored in ascending order.
    function [bool] = ascend( obj )
      bool = obj.ascend_ ;
    end
    
    %> @brief Returns if this TFIMTriangleCircuit is stored in descending order.
    function [bool] = descend( obj )
      bool = ~obj.ascend_ ;
    end
    
    %> @brief Reorders a descending TFIMTriangleCircuit to an ascending
    %> TFIMTriangleCircuit. Involves a copy of all gates. No effect if already
    %> ascending.
    function makeAscend( obj )
      if ~obj.ascend_
        gatesCopy = copy( obj.gates_ );
        for Di = 1:obj.nbGates % loop over all gates
          [ld, id] = obj.descTo2DIndex(obj.nbQubits, Di);
          [Ai] = obj.descAscLinearIndex(obj.nbQubits, ld, id);
          obj.gates_( Ai ) = gatesCopy( Di );
        end
        obj.ascend_ = true;
      end
    end
    
    %> @brief Reorders an ascending TFIMTriangleCircuit to descending
    %> TFIMTriangleCircuit. Involves a copy of all gates. No effect if already
    %> descending.
    function makeDescend( obj )
      if obj.ascend_
        gatesCopy = copy( obj.gates_ );
        for Ai = 1:obj.nbGates % loop over all gates
          [la, ia] = obj.ascTo2DIndex(obj.nbQubits, Ai);
          [Di] = obj.ascDescLinearIndex(obj.nbQubits, la, ia);
          obj.gates_( Di ) = gatesCopy( Ai );
        end
        obj.ascend_ = false;
      end
    end
    
    function merge( obj, gate, side )
      % verify if gate can be merged with triangle
      nbQubits = double(obj.nbQubits);
      if isa(gate,'qclab.qgates.QRotationGate1')
        assert( isa( gate, class(obj.gates_( 1 ) ) ) ); % single-spin gate
        q = double(gate.qubit);
        assert( q < nbQubits );
        ia = 2*(nbQubits - 1 - q);
        id = 2*q;
      elseif isa(gate,'qclab.qgates.QRotationGate2')
        assert( isa( gate, class(obj.gates_( 2 ) ) ) ); % two-spin gate
        q = double(gate.qubit);
        assert( q < nbQubits-1 );
        ia = 2*(nbQubits - 1 - q)-1;
        id = 2*q + 1;
      else
        error('incompatible gate type');
      end
        
      if obj.ascend % ascending triangle
        if strcmp( side, 'L' ) % gate on left side
          currentLayer = 0;
          for i = ia-1:-1:0
            linearIndex = obj.ascAscLinearIndex( nbQubits, currentLayer, i );
            [ obj.gates_(linearIndex), obj.gates_(linearIndex+1), gate ] = ...
             turnover( gate, obj.gates_(linearIndex), obj.gates_(linearIndex+1) );
            currentLayer = currentLayer + 1;
          end
          % fuse from left
          linearIndex = obj.ascAscLinearIndex( nbQubits, currentLayer, 0 );
          obj.gates_(linearIndex).fuse( gate, 'L' );
        else % gate on right side
          currentLayer = 2*nbQubits-2;
          for i = id:2*nbQubits-3
            linearIndex1 = obj.descAscLinearIndex( nbQubits, currentLayer, id );
            linearIndex2 = obj.descAscLinearIndex( nbQubits, currentLayer, id + 1);
            [ gate, obj.gates_(linearIndex1), obj.gates_(linearIndex2) ] = ...
             turnover( obj.gates_(linearIndex1), obj.gates_(linearIndex2), gate );
            currentLayer = currentLayer - 1;
          end
          % fuse from right side
          linearIndex = obj.descAscLinearIndex( nbQubits, currentLayer, id );
          obj.gates_(linearIndex).fuse( gate, 'R' );
        end
      else % descending triangle
        if strcmp( side, 'L' ) % gate on left side
          currentLayer = 0;
          for i = ia-1:-1:0
            linearIndex1 = obj.ascDescLinearIndex( nbQubits, currentLayer, i );
            linearIndex2 = obj.ascDescLinearIndex( nbQubits, currentLayer, i + 1 );
            [ obj.gates_(linearIndex1), obj.gates_(linearIndex2), gate ] = ...
             turnover( gate, obj.gates_(linearIndex1), obj.gates_(linearIndex2) );
            currentLayer = currentLayer + 1;
          end
          % fuse from the left
          linearIndex = obj.ascDescLinearIndex( nbQubits, currentLayer, 0 );
          obj.gates_(linearIndex).fuse( gate, 'L' );
        else % gate on right side
          currentLayer = 2*nbQubits-2;
          for i = id:2*nbQubits-3
            linearIndex = obj.descDescLinearIndex( nbQubits, currentLayer, id );
            [ gate, obj.gates_(linearIndex), obj.gates_(linearIndex + 1) ] = ...
             turnover( obj.gates_(linearIndex), obj.gates_(linearIndex + 1), gate );
            currentLayer = currentLayer - 1;
          end
          % fuse from right side
          linearIndex = obj.descDescLinearIndex( nbQubits, currentLayer, id );
          obj.gates_(linearIndex).fuse( gate, 'R' );
        end
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
    
     %> @brief Conversion to a TFIMSquareCircuit `TFIMSqC`. Doesn't modify `obj`.
    function TFIMSqC = f3c.TFIMSquareCircuit( obj )
      nbQubits = double( obj.nbQubits );
      C = qclab.QCircuit( nbQubits );
      gatesNew = repmat( qclab.qgates.Identity, 1, obj.nbGates );
      
      if obj.descend
        obj.makeAscend;
      end
      
      gatesCopy = copy( obj.gates_ );
      
      %Turnover the even slanted layers from right to left through Triangle
      for turnoverLayer = 2*nbQubits-2:-2:2
        for ia_init = 0:2*nbQubits-2-turnoverLayer
          ia = ia_init;
          % get the gate from the Triangle
          [linearIndex] = obj.ascAscLinearIndex( nbQubits, turnoverLayer, ia );
          gate = gatesCopy( linearIndex );
          for l = turnoverLayer-1:-1:0
            linearIndex = obj.ascAscLinearIndex( nbQubits, l, ia );
            [gate, gatesCopy(linearIndex), gatesCopy(linearIndex+1)] = ...
             turnover(gatesCopy(linearIndex), gatesCopy(linearIndex+1), gate);
            ia = ia + 1;
          end
          coordSum1 = nbQubits - 1 - turnoverLayer/2;
          coordSum2 = coordSum1 + nbQubits - 1;
          if mod(ia_init,2) == 0 % single-spin gate
            SqIndex = f3c.TFIMSquareCircuit.linearIndex( nbQubits, ia_init/2, ...
                                                      coordSum1 - (ia_init/2) );
          else % two-spin gate
            SqIndex = f3c.TFIMSquareCircuit.linearIndex( nbQubits, ...
                                   (ia_init-1)/2, coordSum2 - ((ia_init-1)/2) );
          end
          gatesNew( SqIndex ) = gate;
        end
      end
      
      % Add the remaining layers without moving
      coordSum1 = nbQubits - 1;
      coordSum2 = coordSum1 + nbQubits - 1;
      for ia=0:2*nbQubits-2
        linearIndex = obj.ascAscLinearIndex(nbQubits, 0, ia);
        if mod(ia,2) == 0 % single-spin gate
          SqIndex = f3c.TFIMSquareCircuit.linearIndex( nbQubits, ia/2, ...
                                                           coordSum1 - (ia/2) );
        else % two-spin gate
          SqIndex = f3c.TFIMSquareCircuit.linearIndex( nbQubits, (ia-1)/2, ...
                                                       coordSum2 - ((ia-1)/2) );
        end
        gatesNew( SqIndex ) = gatesCopy( linearIndex );
      end
      
      for moveLayer = 1:2:2*nbQubits-3
        coordSum1 = (moveLayer-1)/2 + nbQubits;
        coordSum2 = coordSum1 + nbQubits - 1;
        for ia=0:2*nbQubits-2-moveLayer
          linearIndex = obj.ascAscLinearIndex(nbQubits, moveLayer, ia);
          if mod(ia,2) == 0 % single-spin gate
            l = (ia + moveLayer + 1)/2;
            SqIndex = f3c.TFIMSquareCircuit.linearIndex( nbQubits, l, ...
                                                                coordSum1 - l );
          else % two-spin gate
            l = (ia + moveLayer)/2;
            SqIndex = f3c.TFIMSquareCircuit.linearIndex( nbQubits, l, ...
                                                                coordSum2 - l );
          end
          gatesNew( SqIndex ) = gatesCopy( linearIndex );
        end
      end
      
      % assign output
      C.insert( 1:length(gatesNew), gatesNew );
      TFIMSqC = f3c.TFIMSquareCircuit( C );
    end
  end
  
  methods (Static)
    
    % Disable QCircuit functionality to limit modifiability of 
    % TFIMTriangleCircuits
    function setOffset(obj, offset)
      assert( false );
    end
    
    function [gates] = gateHandles( obj, pos )
      assert( false );
    end
    
    %
    % Override QCircuit Modifiers for SquareCircuit
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
    
    % Descending linear index from Descending 2D indices
    function [Di] = descDescLinearIndex(nbQubits, ld, id)
      assert( ld < 2*nbQubits - 1);
      assert( id <= ld );
      Di = (ld * (ld+1)) / 2 + id + 1;  
    end
    
    % Descending linear index from Ascending 2D indices
    function [Di] = ascDescLinearIndex(nbQubits, la, ia)
      assert( la < 2*nbQubits - 1 );
      assert( ia < 2*nbQubits - 1 - la );
      Di = ((la+1)*(la+2) + ia*(ia + 2*la + 1))/2;
    end
    
    % Ascending linear index from Ascending 2D indices
    function [Ai] = ascAscLinearIndex(nbQubits, la, ia)
      assert( la < 2*nbQubits - 1 );
      assert( ia < 2*nbQubits - 1 - la );
      nbQubits = double(nbQubits);
      Ai = ia + 1 + la * (2*nbQubits - 1 - (la - 1)/2); 
    end
    
    % Ascending linear index from Descending 2D indices
    function [Ai] = descAscLinearIndex(nbQubits, ld, id)
      assert( ld < 2*nbQubits - 1);
      assert( id <= ld );
      nbQubits = double(nbQubits);
      Ai = ld + 1 + ((4*nbQubits - 3 - id)*id)/2;
    end
    
    % Compute 2D descending indices from descending linear index
    function [ld, id] = descTo2DIndex(nbQubits, Di)
      assert( Di <= nbQubits * (2*nbQubits - 1) );
      % descending layer index
      ld = floor( (sqrt(8*Di - 7) - 1)/2 );
      % descending index within layer
      id = Di - 1 - (ld*(ld+1))/2;
    end
    
    % Compute 2D ascending indices from ascending linear index
    function [la, ia] = ascTo2DIndex(nbQubits, Ai)
      assert( Ai <= nbQubits * (2*nbQubits - 1) );
      nbQubits = double(nbQubits);
      % ascending layer index
      la = floor( (4 * nbQubits - 1 - ...
                       sqrt( 8 *( nbQubits * (2*nbQubits-1) - Ai) + 9 ) ) / 2 );
      % ascending index within layer
      ia = Ai - 1 - la * (2 * nbQubits - 1 - (la - 1)/2);
    end
  end
end