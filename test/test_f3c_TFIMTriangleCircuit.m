classdef test_f3c_TFIMTriangleCircuit < matlab.unittest.TestCase
  methods (Test)
    function test_TFIMTriangleCircuit_index_functions( test)
       % Create an ascending TFIM Triangle Circuit on 4 spins.
      AC = test.randomTFIMAscTr(4);      
      ascend = true;
      ATrC = f3c.TFIMTriangleCircuit( AC, ascend );
      
      % qubits
      test.verifyEqual( ATrC.nbQubits, int32(4) );
      test.verifyEqual( ATrC.qubit, int32(0) );
      test.verifyEqual( ATrC.qubits, int32(0:3) );
      
      % nbGates
      test.verifyEqual( ATrC.nbGates, 28 );
      
      % offset
      test.verifyEqual( ATrC.offset, int32(0) );
      
      % ascend / descend
      test.verifyTrue( ATrC.ascend );
      test.verifyFalse( ATrC.descend );
      
      % matrix
      test.verifyEqual( ATrC.matrix, AC.matrix, 'AbsTol', eps(100) );
      
      % Test ascend index function
      for i = 1:ATrC.nbGates
        [la, ia] = ATrC.ascTo2DIndex( ATrC.nbQubits,  i );
        test.verifyEqual( ATrC.ascAscLinearIndex( ATrC.nbQubits, la, ia ), i );
      end
      
      % Reorder to Descending
      Amat = ATrC.matrix;
      ATrC.makeDescend;
      test.verifyFalse( ATrC.ascend );
      test.verifyTrue( ATrC.descend );
      test.verifyEqual( ATrC.matrix, Amat, 'AbsTol', eps(100) );
      
      % Create a descending TFIM Triangle Circuit on 4 spins.
      DC = test.randomTFIMDescTr(4);      
      ascend = false;
      DTrC = f3c.TFIMTriangleCircuit( DC, ascend );
      
      % qubits
      test.verifyEqual( DTrC.nbQubits, int32(4) );
      test.verifyEqual( DTrC.qubit, int32(0) );
      test.verifyEqual( DTrC.qubits, int32(0:3) );
      
      % nbGates
      test.verifyEqual( DTrC.nbGates, 28 );
      
      % offset
      test.verifyEqual( DTrC.offset, int32(0) );
      
      % ascend / descend
      test.verifyFalse( DTrC.ascend );
      test.verifyTrue( DTrC.descend );
      
      % matrix
      test.verifyEqual( DTrC.matrix, DC.matrix, 'AbsTol', eps(100) );
      
      % Test descend index functions
      for i = 1:DTrC.nbGates
        [ld, id] = DTrC.descTo2DIndex( DTrC.nbQubits,  i );
        test.verifyEqual( DTrC.descDescLinearIndex( DTrC.nbQubits, ld, id ), i );
      end
      
      % Reorder to Ascending
      Dmat = DTrC.matrix;
      DTrC.makeAscend;
      test.verifyTrue( DTrC.ascend );
      test.verifyFalse( DTrC.descend );
      test.verifyEqual( DTrC.matrix, Dmat, 'AbsTol', eps(100) );
      
    end
    
    function test_TFIMTriangleCircuit_toSquare( test )
      % Transform Ascending Triangle to Square
      AC = test.randomTFIMAscTr(5);      
      ascend = true;
      ATrC = f3c.TFIMTriangleCircuit( AC, ascend );
      ATrCmat = ATrC.matrix;
      SqC = f3c.TFIMSquareCircuit( ATrC );
      test.verifyEqual( SqC.matrix, ATrCmat, 'AbsTol', eps(1000) );
      
      % Transform Descending Triangle to Square
      DC = test.randomTFIMDescTr(6);
      ascend = false;
      DTrC = f3c.TFIMTriangleCircuit( DC, ascend );
      DTrCmat = DTrC.matrix;
      SqC = f3c.TFIMSquareCircuit( DTrC );
      test.verifyEqual( SqC.matrix, DTrCmat, 'AbsTol', eps(1000) );
    end
    
    function test_TFIMTriangleCircuit_mergeGate( test )
      Z = @qclab.qgates.RotationZ;
      XX = @qclab.qgates.RotationXX;
      
      % Ascending Triangle, gates on left hand side
      AC = test.randomTFIMAscTr(5);      
      ascend = true;
      ATrC = f3c.TFIMTriangleCircuit( AC, ascend );
      ATrCmat = ATrC.matrix;
      ATrCO = copy( ATrC );
      
      % Single-spin gates
      % q0
      gate = Z( 0, randn );
      gatemat = kron( gate.matrix, qclab.qId(4) );
      ATrC.merge( gate, 'L' );
      test.verifyEqual( ATrC.matrix, ATrCmat * gatemat, 'AbsTol', eps(1000) );
      
      % q1
      ATrC = copy( ATrCO );
      gate = Z( 1, randn );
      gatemat = kron( kron(qclab.qId(1), gate.matrix) , qclab.qId(3) );
      ATrC.merge( gate, 'L' );
      test.verifyEqual( ATrC.matrix, ATrCmat * gatemat, 'AbsTol', eps(1000) );
      
      % q2
      ATrC = copy( ATrCO );
      gate = Z( 2, randn );
      gatemat = kron( kron(qclab.qId(2), gate.matrix) , qclab.qId(2) );
      ATrC.merge( gate, 'L' );
      test.verifyEqual( ATrC.matrix, ATrCmat * gatemat, 'AbsTol', eps(1000) );
      
      % q3
      ATrC = copy( ATrCO );
      gate = Z( 3, randn );
      gatemat = kron( kron(qclab.qId(3), gate.matrix) , qclab.qId(1) );
      ATrC.merge( gate, 'L' );
      test.verifyEqual( ATrC.matrix, ATrCmat * gatemat, 'AbsTol', eps(1000) );
      
      % q4
      ATrC = copy( ATrCO );
      gate = Z( 4, randn );
      gatemat = kron(qclab.qId(4), gate.matrix);
      ATrC.merge( gate, 'L' );
      test.verifyEqual( ATrC.matrix, ATrCmat * gatemat, 'AbsTol', eps(1000) );
      
      % Two-spin gates
      % q0,q1
      ATrC = copy( ATrCO );
      gate = XX( [0,1], randn );
      gatemat = kron( gate.matrix, qclab.qId(3) );
      ATrC.merge( gate, 'L' );
      test.verifyEqual( ATrC.matrix, ATrCmat * gatemat, 'AbsTol', eps(1000) );
      
      % q1,q2
      ATrC = copy( ATrCO );
      gate = XX( [1,2], randn );
      gatemat = kron( kron(qclab.qId(1), gate.matrix) , qclab.qId(2) );
      ATrC.merge( gate, 'L' );
      test.verifyEqual( ATrC.matrix, ATrCmat * gatemat, 'AbsTol', eps(1000) );
      
      % q2,q3
      ATrC = copy( ATrCO );
      gate = XX( [2,3], randn );
      gatemat = kron( kron(qclab.qId(2), gate.matrix) , qclab.qId(1) );
      ATrC.merge( gate, 'L' );
      test.verifyEqual( ATrC.matrix, ATrCmat * gatemat, 'AbsTol', eps(1000) );
      
      % q3,q4
      ATrC = copy( ATrCO );
      gate = XX( [3,4], randn );
      gatemat = kron(qclab.qId(3), gate.matrix);
      ATrC.merge( gate, 'L' );
      test.verifyEqual( ATrC.matrix, ATrCmat * gatemat, 'AbsTol', eps(1000) );
      
      % Ascending Triangle, gates on right hand side
      % Single-spin gates
      % q0
      ATrC = copy( ATrCO );
      gate = Z( 0, randn );
      gatemat = kron( gate.matrix, qclab.qId(4) );
      ATrC.merge( gate, 'R' );
      test.verifyEqual( ATrC.matrix, gatemat * ATrCmat, 'AbsTol', eps(1000) );
      
      % q1
      ATrC = copy( ATrCO );
      gate = Z( 1, randn );
      gatemat = kron( kron(qclab.qId(1), gate.matrix) , qclab.qId(3) );
      ATrC.merge( gate, 'R' );
      test.verifyEqual( ATrC.matrix, gatemat * ATrCmat, 'AbsTol', eps(1000) );
      
      % q2
      ATrC = copy( ATrCO );
      gate = Z( 2, randn );
      gatemat = kron( kron(qclab.qId(2), gate.matrix) , qclab.qId(2) );
      ATrC.merge( gate, 'R' );
      test.verifyEqual( ATrC.matrix, gatemat * ATrCmat, 'AbsTol', eps(1000) );
      
      % q3
      ATrC = copy( ATrCO );
      gate = Z( 3, randn );
      gatemat = kron( kron(qclab.qId(3), gate.matrix) , qclab.qId(1) );
      ATrC.merge( gate, 'R' );
      test.verifyEqual( ATrC.matrix, gatemat * ATrCmat, 'AbsTol', eps(1000) );
      
      % q4
      ATrC = copy( ATrCO );
      gate = Z( 4, randn );
      gatemat = kron(qclab.qId(4), gate.matrix);
      ATrC.merge( gate, 'R' );
      test.verifyEqual( ATrC.matrix, gatemat * ATrCmat, 'AbsTol', eps(1000) );
      
      % Two-spin gates
      % q0,q1
      ATrC = copy( ATrCO );
      gate = XX( [0,1], randn );
      gatemat = kron( gate.matrix, qclab.qId(3) );
      ATrC.merge( gate, 'R' );
      test.verifyEqual( ATrC.matrix, gatemat * ATrCmat, 'AbsTol', eps(1000) );
      
      % q1,q2
      ATrC = copy( ATrCO );
      gate = XX( [1,2], randn );
      gatemat = kron( kron(qclab.qId(1), gate.matrix) , qclab.qId(2) );
      ATrC.merge( gate, 'R' );
      test.verifyEqual( ATrC.matrix, gatemat * ATrCmat, 'AbsTol', eps(1000) );
      
      % q2,q3
      ATrC = copy( ATrCO );
      gate = XX( [2,3], randn );
      gatemat = kron( kron(qclab.qId(2), gate.matrix) , qclab.qId(1) );
      ATrC.merge( gate, 'R' );
      test.verifyEqual( ATrC.matrix, gatemat * ATrCmat , 'AbsTol', eps(1000) );
      
      
      % Descending triangle, gates on left hand side
      DC = test.randomTFIMDescTr(5);
      ascend = false;
      DTrC = f3c.TFIMTriangleCircuit( DC, ascend );
      DTrCmat = DTrC.matrix;
      DTrCO = copy( DTrC );
      
      % Single-spin gates
      % q0
      gate = Z( 0, randn );
      gatemat = kron( gate.matrix, qclab.qId(4) );
      DTrC.merge( gate, 'L' );
      test.verifyEqual( DTrC.matrix, DTrCmat * gatemat, 'AbsTol', eps(1000) );
      
      % q1
      DTrC = copy( DTrCO );
      gate = Z( 1, randn );
      gatemat = kron( kron(qclab.qId(1), gate.matrix) , qclab.qId(3) );
      DTrC.merge( gate, 'L' );
      test.verifyEqual( DTrC.matrix, DTrCmat * gatemat, 'AbsTol', eps(1000) );
      
      % q2
      DTrC = copy( DTrCO );
      gate = Z( 2, randn );
      gatemat = kron( kron(qclab.qId(2), gate.matrix) , qclab.qId(2) );
      DTrC.merge( gate, 'L' );
      test.verifyEqual( DTrC.matrix, DTrCmat * gatemat, 'AbsTol', eps(1000) );
      
      % q3
      DTrC = copy( DTrCO );
      gate = Z( 3, randn );
      gatemat = kron( kron(qclab.qId(3), gate.matrix) , qclab.qId(1) );
      DTrC.merge( gate, 'L' );
      test.verifyEqual( DTrC.matrix, DTrCmat * gatemat, 'AbsTol', eps(1000) );
      
      % q4
      DTrC = copy( DTrCO );
      gate = Z( 4, randn );
      gatemat = kron(qclab.qId(4), gate.matrix);
      DTrC.merge( gate, 'L' );
      test.verifyEqual( DTrC.matrix, DTrCmat * gatemat, 'AbsTol', eps(1000) );
      
      % Two-spin gates
      % q0,q1
      DTrC = copy( DTrCO );
      gate = XX( [0,1], randn );
      gatemat = kron( gate.matrix, qclab.qId(3) );
      DTrC.merge( gate, 'L' );
      test.verifyEqual( DTrC.matrix, DTrCmat * gatemat, 'AbsTol', eps(1000) );
      
      % q1,q2
      DTrC = copy( DTrCO );
      gate = XX( [1,2], randn );
      gatemat = kron( kron(qclab.qId(1), gate.matrix) , qclab.qId(2) );
      DTrC.merge( gate, 'L' );
      test.verifyEqual( DTrC.matrix, DTrCmat * gatemat, 'AbsTol', eps(1000) );
      
      % q2,q3
      DTrC = copy( DTrCO );
      gate = XX( [2,3], randn );
      gatemat = kron( kron(qclab.qId(2), gate.matrix) , qclab.qId(1) );
      DTrC.merge( gate, 'L' );
      test.verifyEqual( DTrC.matrix, DTrCmat * gatemat, 'AbsTol', eps(1000) );
      
      % q3,q4
      DTrC = copy( DTrCO );
      gate = XX( [3,4], randn );
      gatemat = kron(qclab.qId(3), gate.matrix);
      DTrC.merge( gate, 'L' );
      test.verifyEqual( DTrC.matrix, DTrCmat * gatemat, 'AbsTol', eps(1000) );
      
      % Descending Triangle, gates on right hand side
      % Single-spin gates
      % q0
      DTrC = copy( DTrCO );
      gate = Z( 0, randn );
      gatemat = kron( gate.matrix, qclab.qId(4) );
      DTrC.merge( gate, 'R' );
      test.verifyEqual( DTrC.matrix, gatemat * DTrCmat, 'AbsTol', eps(1000) );
      
      % q1
      DTrC = copy( DTrCO );
      gate = Z( 1, randn );
      gatemat = kron( kron(qclab.qId(1), gate.matrix) , qclab.qId(3) );
      DTrC.merge( gate, 'R' );
      test.verifyEqual( DTrC.matrix, gatemat * DTrCmat, 'AbsTol', eps(1000) );
      
      % q2
      DTrC = copy( DTrCO );
      gate = Z( 2, randn );
      gatemat = kron( kron(qclab.qId(2), gate.matrix) , qclab.qId(2) );
      DTrC.merge( gate, 'R' );
      test.verifyEqual( DTrC.matrix, gatemat * DTrCmat, 'AbsTol', eps(1000) );
      
      % q3
      DTrC = copy( DTrCO );
      gate = Z( 3, randn );
      gatemat = kron( kron(qclab.qId(3), gate.matrix) , qclab.qId(1) );
      DTrC.merge( gate, 'R' );
      test.verifyEqual( DTrC.matrix, gatemat * DTrCmat, 'AbsTol', eps(1000) );
      
      % q4
      DTrC = copy( DTrCO );
      gate = Z( 4, randn );
      gatemat = kron(qclab.qId(4), gate.matrix);
      DTrC.merge( gate, 'R' );
      test.verifyEqual( DTrC.matrix, gatemat * DTrCmat, 'AbsTol', eps(1000) );
      
      % Two-spin gates
      % q0,q1
      DTrC = copy( DTrCO );
      gate = XX( [0,1], randn );
      gatemat = kron( gate.matrix, qclab.qId(3) );
      DTrC.merge( gate, 'R' );
      test.verifyEqual( DTrC.matrix, gatemat * DTrCmat, 'AbsTol', eps(1000) );
      
      % q1,q2
      DTrC = copy( DTrCO );
      gate = XX( [1,2], randn );
      gatemat = kron( kron(qclab.qId(1), gate.matrix) , qclab.qId(2) );
      DTrC.merge( gate, 'R' );
      test.verifyEqual( DTrC.matrix, gatemat * DTrCmat, 'AbsTol', eps(1000) );
      
      % q2,q3
      DTrC = copy( DTrCO );
      gate = XX( [2,3], randn );
      gatemat = kron( kron(qclab.qId(2), gate.matrix) , qclab.qId(1) );
      DTrC.merge( gate, 'R' );
      test.verifyEqual( DTrC.matrix, gatemat * DTrCmat , 'AbsTol', eps(1000) );
      
      % q3,q4
      DTrC = copy( DTrCO );
      gate = XX( [3,4], randn );
      gatemat = kron(qclab.qId(3), gate.matrix);
      DTrC.merge( gate, 'R' );
      test.verifyEqual( DTrC.matrix, gatemat * DTrCmat, 'AbsTol', eps(1000) );
      
    end
  end
  
  methods (Static)
    function C = randomTFIMAscTr( n )
      Z = @qclab.qgates.RotationZ;
      XX = @qclab.qgates.RotationXX;
      
      C = qclab.QCircuit( n );
      
      for l = 0:2*n-2
        q = n-1;
        for ia = 0:2*n-2-l
          if mod(ia,2) == 0 % single-spin
            C.push_back( Z( q, randn ) );
            q = q - 1;
          else % two-spin
            C.push_back( XX( [q, q+1], randn ) );
          end
        end
      end
      
    end
    
    function C = randomTFIMDescTr( n )
      Z = @qclab.qgates.RotationZ;
      XX = @qclab.qgates.RotationXX;
      
      C = qclab.QCircuit( n );
      
      for l = 0:2*n-2
        if mod(l,2) == 0 % even layer
          q = n - 1 - l/2;
          for id = 0:l
            if mod(id,2) == 0 % even id
              C.push_back( Z( q, randn ) );
            else % odd id
              C.push_back( XX( [q,q+1], randn ) );
              q = q + 1;
            end
          end
        else % odd layer
          q = n - 1 - (l+1)/2;
          for id = 0:l
            if mod(id,2) == 0 % even id
              C.push_back( XX( [q,q+1], randn ) );
              q = q + 1;
            else % odd id
              C.push_back( Z( q, randn ) );
            end
          end
        end
      end
    end
  end
end