%> @file TFIMSquareCircuit.m
%> @brief Implements TFIMSquareCircuit class.
% ==============================================================================
%> @class TFIMSquareCircuit
%> @brief Class for representing a square quantum circuit for Hamiltonian
%> simulation of TFIM models
%>
%> In a TFIMSquareCircuit, the one- and two-spin gates appear in the following
%> pattern:
%>
%>    O |     O |     O |     O |         01  05  08  12  15  19  22  26
%>    O | |   O | |   O | |   O | |       02  06  09  13  16  20  23  27
%>    O   | | O   | | O   | | O   | |     03  07  10  14  17  21  24  28
%>    O     | O     | O     | O     |     04      11      18      25
%>
%> where:
%>
%>    | : two-spin rotation gate
%>    |
%>
%>    O : one-spin rotation gate (transverse-field) 
%>
%>
%>  A single layer is given by:
%>
%>  O |           01  05
%>  O | |         02  06
%>  O   | |       03  07
%>  O     |       04
%>
%>
%> An N-qubit TFIMSquareCircuit consists of N layers with 2N-1 gates each (N
%> single-spin and N-1 two-spin)
%
% (C) Copyright Daan Camps and Roel Van Beeumen 2021.  
% ==============================================================================
classdef TFIMSquareCircuit < qclab.QCircuit
  
  methods
    % Constructor ==============================================================
    %> Constructor for TFIMSquareCircuit objects
    %>
    %>   1. TFIMSquareCircuit() : create a 1-qubit, empty SquareCircuit
    %>
    %>   2. TFIMSquareCircuit( circuit ) : create a TFIMSquareCircuit from a
    %>      regular circuit of correct dimensions.
    %
    % ==========================================================================
    function obj = TFIMSquareCircuit( circuit )
      obj = obj@qclab.QCircuit( 1, 0 );
      if nargin == 0
        obj.push_back( qclab.qgates.RotationZ(0, 0 ) );
      else % initialize with circuit
        obj.update( circuit );
      end
    end
    
    %> @brief Update this TFIMSquareCircuit with a QCircuit of right dimensions.
    function update( obj, circuit )
      assert( isa( circuit, 'qclab.QCircuit' ) );
      nbQubits = double( circuit.nbQubits );
      nbGates = circuit.nbGates;
      assert( nbQubits * (2*nbQubits - 1) == nbGates );
      type1 = class( circuit.gateHandles(1) ); % 1-spin rotation
      type2 = class( circuit.gateHandles( nbQubits + 1 ) ); % 2-spin rotation
      i = 1;
      for l = 1:nbQubits
        % one-spin rotations
        for q = 0:nbQubits-1
          gatei = circuit.gateHandles(i);
          i = i + 1;
          assert( isa( gatei, type1 ) );
          assert( gatei.qubit == q );
        end
        % two-spin rotations
        for q = 0:nbQubits-2
          gatei = circuit.gateHandles(i);
          i = i + 1;
          assert( isa( gatei, type2 ) );
          assert( gatei.qubit == q );
        end
      end
      % all OK - initialize gates
      obj.nbQubits_ = nbQubits ;
      obj.gates_ = circuit.gateHandles() ;
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
    
    %> @brief Conversion to an ascending TFIMTriangleCircuit `TFIMTrC`
    function TFIMTrC = f3c.TFIMTriangleCircuit( obj )
      nbQubits = double( obj.nbQubits );
      C = qclab.QCircuit( nbQubits );
      gatesNew = repmat( qclab.qgates.Identity, 1, obj.nbGates );
      gatesCopy = copy( obj.gates_ );

      % First fixed layer
      coordSum1 = nbQubits - 1;
      for l = 0:nbQubits-2 % loop over vertical layers in square
        linearIndex = obj.linearIndex( nbQubits, l, coordSum1 - l );
        TrIndex = f3c.TFIMTriangleCircuit.ascAscLinearIndex( nbQubits, 0, 2*l );
        gatesNew( TrIndex ) = gatesCopy( linearIndex );
        gatesNew( TrIndex + 1 ) = gatesCopy( linearIndex + nbQubits - 1 );
      end
      linearIndex = obj.linearIndex( nbQubits, nbQubits - 1, ...
                                                      coordSum1 - nbQubits + 1);
      gatesNew( TrIndex + 2 ) = gatesCopy( linearIndex );
      
      % Other move layers
      for moveLayer = 1:2:2*nbQubits - 2 % loop over layers in triangle
         coordSum1 = nbQubits - 1 + (moveLayer+1)/2;
         j = 0;
         for l = (moveLayer+1)/2:nbQubits-1
           linearIndex = obj.linearIndex( nbQubits, l, coordSum1 - l );
           TrIndex = f3c.TFIMTriangleCircuit.ascAscLinearIndex( nbQubits, ...
                                                                  moveLayer, j);
           j = j + 1;
           gatesNew( TrIndex ) = gatesCopy( linearIndex );
           if coordSum1 - l > 0
            gatesNew( TrIndex + 1 ) = gatesCopy( linearIndex + nbQubits - 1 );
            j = j + 1;
           end
         end
       end
      
      % Turnover layers
      for turnoverLayer = nbQubits-2:-1:0 % loop over top half of square
        ia_init = 2*nbQubits - 3;
        k = 0;  
        for l = turnoverLayer:-1:0 % loop within layer         
          % first (single-spin) gate
          ia = ia_init - k;      
          linearIndex = obj.linearIndex( nbQubits, l, turnoverLayer - l );
          gate = gatesCopy( linearIndex );          
          for i = 0:2*(nbQubits - turnoverLayer - 1)-1 % over triangle layers
            linearIndex = f3c.TFIMTriangleCircuit.ascAscLinearIndex( nbQubits, ...
                                                                        i, ia );
            [ gatesNew(linearIndex), gatesNew(linearIndex + 1), gate] = ...
             turnover( gate, gatesNew(linearIndex), gatesNew(linearIndex + 1) );
            ia = ia - 1;
          end
          ia = ia + 1;
          linearIndex = f3c.TFIMTriangleCircuit.ascAscLinearIndex( nbQubits, ...
                                           2*(nbQubits - turnoverLayer-1), ia ); 
          gatesNew( linearIndex ) = gate;
          k = k + 1;
          
          % second (two-spin) gate
          if l > 0
            ia = ia_init - k; 
            linearIndex = obj.linearIndex( nbQubits, l - 1, ...
                                                 turnoverLayer - l + nbQubits );
            gate = gatesCopy( linearIndex );
            for i = 0:2*(nbQubits - turnoverLayer - 1)-1 % over triangle layers
              linearIndex = f3c.TFIMTriangleCircuit.ascAscLinearIndex( nbQubits, ...
                                                                        i, ia );
              [ gatesNew(linearIndex), gatesNew(linearIndex + 1), gate] = ...
                turnover( gate, gatesNew(linearIndex), gatesNew(linearIndex + 1) );
              ia = ia - 1;
            end
            ia = ia + 1;
            linearIndex = f3c.TFIMTriangleCircuit.ascAscLinearIndex( nbQubits,...
                                           2*(nbQubits - turnoverLayer-1), ia );
            gatesNew( linearIndex ) = gate;
            k = k + 1;
          end
        end
      end
      
      % assign output
      C.insert( 1:length(gatesNew), gatesNew );
      TFIMTrC = f3c.TFIMTriangleCircuit( C, true );
    end
  end
  
  methods (Static)
    
    % Disable QCircuit functionality to limit modifiability of TFIMSquareCircuits
    
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
    
    % Compute linear index from 2D indices. 
    % Layer index: 0 <= l < nbQubits
    % j index:      0 <= j < 2*nbQubits-1
    function [idx] = linearIndex(nbQubits, l, j )
      assert( qclab.isNonNegInteger(j) && qclab.isNonNegInteger(l) );
      assert( l < nbQubits );
      assert( j < 2*nbQubits-1 );
      nbQubits = double(nbQubits);
      idx = l * (2 * nbQubits - 1) + j + 1;
    end
    
    % Compute 2D indices from linear index.
    function [l, j] = linearTo2DIndex(nbQubits, idx)
      assert( qclab.isNonNegInteger(idx - 1) );
      assert( idx <= nbQubits * (2*nbQubits - 1) );
      nbQubits = double(nbQubits);
      idx = idx - 1;
      l = floor(idx / (2 * nbQubits - 1));
      j = idx - l * (2 * nbQubits - 1);
    end
  end
end %TFIMSquareCircuit