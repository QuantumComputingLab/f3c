%> @file SquareCircuit.m
%> @brief Implements SquareCircuit class.
% ==============================================================================
%> @class SquareCircuit
%> @brief Class for representing a square quantum circuit for Hamiltonian
%>  simulation
%>
%> In a SquareCircuit, the gates appear in the following pattern and ordering:
%>
%>  A) Even nbQubits
%>      
%>      |   |   |           01      06      11
%>        |   |   |             04      09      14
%>      |   |   |           02      07      12
%>        |   |   |             05      10      15
%>      |   |   |           03      08      13
%>
%>      
%>  B) Odd nbQubits
%>      
%>      |   |   |   |       01      07      13      19
%>        |   |   |             04      10      16
%>      |   |   |   |       02      08      14      20
%>        |   |   |             05      11      17
%>      |   |   |   |       03      09      15      21
%>        |   |   |             06      12      18
%>
%>
% (C) Copyright Daan Camps and Roel Van Beeumen 2021.  
% ==============================================================================
classdef SquareCircuit < qclab.QCircuit
  
  methods
    
    % Constructor ==============================================================
    %> Constructor for SquareCircuit objects
    %>
    %>   1. SquareCircuit() : create a 1-qubit, empty SquareCircuit
    %>
    %>   2. SquareCircuit( circuit ) : create a SquareCircuit from a regular
    %>      circuit of correct dimensions.
    %
    % ==========================================================================
    function obj = SquareCircuit( circuit )
      obj = obj@qclab.QCircuit( 1, 0 );
      if nargin == 1 % initialize with circuit
          obj.update( circuit );
      end
    end
    
    %> @brief Update this SquareCircuit with a QCircuit of right dimensions.
    function update( obj, circuit )
      assert( isa( circuit, 'qclab.QCircuit' ) );
      nbQubits = double( circuit.nbQubits );
      nbGates = circuit.nbGates ;
      assert( (nbQubits * (nbQubits - 1))/2 == nbGates );
      % assert that all the gates are ordered as expected and of compatible
      % types.
      types = cell( nbQubits - 1, 1 );
      for i = 1:nbQubits - 1
        [qExpected,~] = obj.computeQLIndex(nbQubits, i);
        gatei = circuit.gateHandles(i);
        assert( isa( gatei, 'qclab.qgates.QGate2' ) );
        qubitsi = gatei.qubits ;
        assert( qubitsi(1) == qExpected );
        assert( qubitsi(1) + 1 == qubitsi(2) );
        types{qExpected + 1} = class(gatei);
        
      end
      for i = nbQubits:nbGates
        [qExpected,~] = obj.computeQLIndex(nbQubits, i);
        gatei = circuit.gateHandles(i) ;
        qubitsi = gatei.qubits ;
        assert( qubitsi(1) == qExpected );
        assert( qubitsi(1) + 1 == qubitsi(2) );
        assert( isa( gatei, types{qExpected + 1} ) );
      end
      % all OK - initialize gates
      obj.nbQubits_ = nbQubits ;
      obj.gates_ = circuit.gateHandles() ;
    end
    
    %> @brief Compute the linear index from a qubit and layer index in a 
    %> SquareCircuit.
    %>
    %> @param obj SquareCircuit
    %> @param qubitIndex qubit index of gate ranging from 0 to obj.nbQubits - 2
    %> @param layerIndex layer index of gate ranging from 0 to obj.nbQubits - 1,
    %>                   qubit and layer index must have same even/odd parity.
    %>
    %> @retval linearIndex linear index of gate in quantum circuit.
    function [linearIndex] = linearIndex( obj, qubitIndex, layerIndex )
      
      nbQubits = double( obj.nbQubits );
      linearIndex =  obj.computeLinearIndex( nbQubits, qubitIndex, layerIndex );
    end
    
    %> @brief Compute the qubit and layer index of a gate in a SquareCircuit
    %> from the linear index.
    %>
    %> @param obj SquareCircuit
    %> @param linearIndex linear gate index for SquareCircuit between
    %>                    1:obj.nbGates
    %>
    %> @retval qubitIndex qubit index of gate ranging from 0 to obj.nbQubits - 1
    %> @retval layerIndex layer index of gate ranging from 0 to obj.nbQubits - 1
    %>                    with even layers containing gates with even qubits 
    %>                    and odd layers gates on odd qubits.
    function [qubitIndex, layerIndex] = qlIndex( obj, linearIndex )
      nbQubits = double( obj.nbQubits );
      [qubitIndex, layerIndex] = obj.computeQLIndex( nbQubits, linearIndex );
      
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
    
    %> @brief Conversion to an (ascending) TriangleCircuit `TrC`. Doesn't modify
    %> `obj`.
    function TrC = f3c.TriangleCircuit( obj )
      nbQubits = double( obj.nbQubits );
      C = qclab.QCircuit( nbQubits );
      gatesNew = repmat( qclab.qgates.Identity, 1, obj.nbGates );
      gatesCopy = copy( obj.gates_ );
      
      if mod( nbQubits, 2 ) == 0 % even qubits
        coordSum = nbQubits - 2;
        for qubit = nbQubits-2:-1:0
          linearIndex = obj.linearIndex( qubit, coordSum - qubit );
          TrIndex = f3c.TriangleCircuit.computeAscendLinearIndex( nbQubits, ...
                                                                     qubit, 0 );
          gatesNew( TrIndex ) = gatesCopy( linearIndex );
        end
        turnoverLayerRng = 2:2:nbQubits - 2;
        moveLayerRng = 1:2:nbQubits - 2;
      else % odd qubits
        turnoverLayerRng = 1:2:nbQubits - 2;
        moveLayerRng = 0:2:nbQubits - 2;
      end
      
      % add first part without turnovers
      for moveLayer = moveLayerRng
        coordSum = nbQubits + moveLayer - 1;
        for qubit = nbQubits - 2:-1:moveLayer
         linearIndex = obj.linearIndex( qubit, coordSum - qubit );
         TrIndex = f3c.TriangleCircuit.computeAscendLinearIndex( nbQubits,...
                                                             qubit, moveLayer );
         gatesNew( TrIndex ) = gatesCopy( linearIndex );
        end
      end
      
      % fill-in remaining using turnovers
      for turnoverLayer = turnoverLayerRng
        coordSum = nbQubits - turnoverLayer - 2;
        for qubit = 0:coordSum
          linearIndex = obj.linearIndex( qubit, coordSum - qubit );
          gate = gatesCopy( linearIndex );
          for i = 0:turnoverLayer - 1
            % watch out for double in .qubit
            linearIndex = f3c.TriangleCircuit.computeAscendLinearIndex( ...
                                                  nbQubits, gate.qubit + 1, i );
            [ gatesNew(linearIndex), gatesNew(linearIndex + 1), gate] = ...
             turnover( gate, gatesNew(linearIndex), gatesNew(linearIndex + 1) );
          end
          linearIndex = f3c.TriangleCircuit.computeAscendLinearIndex( ...
                                          nbQubits, gate.qubit, turnoverLayer );
          gatesNew( linearIndex ) = gate;
        end
      end
      
      % assign output
      C.insert( 1:length(gatesNew), gatesNew );
      TrC = f3c.TriangleCircuit( C, true );
      
    end
  end
  
  methods (Static)
    
    % Compute circuit indices
    
    %> @brief Compute the linear index from a qubit and layer index for a 
    %> SquareCircuit.
    %>
    %> @param nbQubits number of qubits in SquareCircuit
    %> @param qubitIndex qubit index of gate ranging from 0 to nbQubits - 2
    %> @param layerIndex layer index of gate ranging from 0 to nbQubits - 1,
    %>                   qubit and layer index must have same even/odd parity.
    %>
    %> @retval linearIndex linear index of gate in quantum circuit.
    function [linearIndex] = computeLinearIndex( nbQubits, qubitIndex, layerIndex )
      assert( mod(qubitIndex, 2) == mod(layerIndex, 2) );
      assert( qubitIndex < nbQubits - 1 && layerIndex < nbQubits );
      layerIndex = double( layerIndex );
      if mod(nbQubits, 2) == 0 % even number of qubits
        if mod( qubitIndex, 2 ) == 0 % even qubit and layer
          linearIndex = (qubitIndex / 2) + 1 + ...
                          (layerIndex / 2) * (nbQubits - 1); 
        else % odd qubit and layer
          linearIndex = ((qubitIndex + nbQubits - 1) / 2 ) + 1 + ...
                          floor( layerIndex /2) * (nbQubits - 1);
        end
      else % odd number of qubits
        if mod( qubitIndex, 2 ) == 0 % even qubit and layer
          linearIndex = (qubitIndex / 2) + 1 + ...
                          (layerIndex / 2) * (nbQubits - 1); 
        else % odd qubit and layer
          linearIndex = ((qubitIndex + nbQubits - 2) / 2 ) + 1 + ...
                          floor( layerIndex /2) * (nbQubits - 1);
        end
      end
    end
    
    %> @brief Compute the qubit and layer index of a gate in a SquareCircuit
    %> from the linear index.
    %>
    %> @param nbQubits number of qubits in SquareCircuit
    %> @param linearIndex linear gate index for SquareCircuit between
    %>                    1:(nbQubits * (nbQubits - 1))/2
    %>
    %> @retval qubitIndex qubit index of gate ranging from 0 to nbQubits - 1
    %> @retval layerIndex layer index of gate ranging from 0 to nbQubits - 1
    %>                    with even layers containing gates with even qubits 
    %>                    and odd layers gates on odd qubits.
    function [qubitIndex, layerIndex] = computeQLIndex( nbQubits, linearIndex )
      assert( linearIndex <= nbQubits * (nbQubits - 1) / 2);
      if mod(nbQubits,2) == 0 % even number of qubits
        qubitIndex = mod( 2 * (linearIndex - 1), nbQubits - 1 );
        layerIndex = 2 * ceil( linearIndex / ( nbQubits - 1 ) ) - 1;
        if mod( qubitIndex, 2 ) == 0
          layerIndex = layerIndex - 1;
        end
      else % odd number of qubits
        wrapped = mod( linearIndex, nbQubits - 1 );
        qubitIndex = mod( 2 * ( wrapped - 1), nbQubits - 1 );
        if wrapped > nbQubits/2 || wrapped == 0
          qubitIndex = qubitIndex + 1;
        end
        layerIndex = ceil( linearIndex / (( nbQubits - 1 )/2) ) - 1;
      end
    end
    
    
    % Disable QCircuit functionality to limit modifiability of SquareCircuits
    
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
    
  end
  
end % SquareCircuit