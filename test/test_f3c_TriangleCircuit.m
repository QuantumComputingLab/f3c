classdef test_f3c_TriangleCircuit < matlab.unittest.TestCase
  methods (Test)
    function test_TriangleCircuit_index_functions( test)
      Circuit = @qclab.QCircuit ;
      TrCircuit = @f3c.TriangleCircuit ;
      XY = @f3c.qgates.RotationXY ;
      
      % Even number of qubits
      % Generate circuit
      C = Circuit( 4 );
      
      % layer 0
      C.push_back( XY([2,3]) );
      C.push_back( XY([1,2]) );      
      C.push_back( XY([0,1]) );
      % layer 1
      C.push_back( XY([2,3]) );
      C.push_back( XY([1,2]) );
      % layer 2
      C.push_back( XY([2,3]) );
      
      ascend = true;
      ATrC = TrCircuit( C, ascend );
      
      % qubits
      test.verifyEqual( ATrC.nbQubits, int32(4) );
      test.verifyEqual( ATrC.qubit, int32(0) );
      test.verifyEqual( ATrC.qubits, int32(0:3) );
      
      % nbGates
      test.verifyEqual( ATrC.nbGates, 6 );
      
      % offset
      test.verifyEqual( ATrC.offset, int32(0) );
      
      % ascend / descend
      test.verifyTrue( ATrC.ascend );
      test.verifyFalse( ATrC.descend );
      
      % matrix
      test.verifyEqual( ATrC.matrix, qclab.qId(4), 'AbsTol', eps(100) );
      
      test.verifyEqual( ATrC.ascendLinearIndex( 2, 0 ) , 1 );
      test.verifyEqual( ATrC.ascendLinearIndex( 1, 0 ) , 2 );
      test.verifyEqual( ATrC.ascendLinearIndex( 0, 0 ) , 3 );
      test.verifyEqual( ATrC.ascendLinearIndex( 2, 1 ) , 4 );
      test.verifyEqual( ATrC.ascendLinearIndex( 1, 1 ) , 5 );
      test.verifyEqual( ATrC.ascendLinearIndex( 2, 2 ) , 6 );
      
      % Ascending Triangle
      C = Circuit( 7 );
      % layer 0
      C.push_back( XY([5,6]) );
      C.push_back( XY([4,5]) );      
      C.push_back( XY([3,4]) );
      C.push_back( XY([2,3]) );
      C.push_back( XY([1,2]) );
      C.push_back( XY([0,1]) );
      
      % layer 1      
      C.push_back( XY([5,6]) );
      C.push_back( XY([4,5]) );      
      C.push_back( XY([3,4]) );
      C.push_back( XY([2,3]) );
      C.push_back( XY([1,2]) );
      
      % layer 2
      C.push_back( XY([5,6]) );
      C.push_back( XY([4,5]) );      
      C.push_back( XY([3,4]) );
      C.push_back( XY([2,3]) );
      
      % layer 3
      C.push_back( XY([5,6]) );
      C.push_back( XY([4,5]) );      
      C.push_back( XY([3,4]) );
      
      % layer 4
      C.push_back( XY([5,6]) );
      C.push_back( XY([4,5]) );      
      
      % layer 5
      C.push_back( XY([5,6]) );
      
      ascend = true;
      ATrC = TrCircuit( C, ascend );
      
      % check index functions
      for i = 1:ATrC.nbGates
        [q, l] = ATrC.QLIndex( i );
        test.verifyEqual( ATrC.ascendLinearIndex( q, l ), i );
      end
      
      % Descending Triangle
      C = Circuit( 7 );
      
      % layer 0
      C.push_back( XY([5,6]) );
      
      % layer 1
      C.push_back( XY([4,5]) );      
      C.push_back( XY([5,6]) );
      
      % layer 2
      C.push_back( XY([3,4]) );
      C.push_back( XY([4,5]) );      
      C.push_back( XY([5,6]) );
      
      % layer 3
      C.push_back( XY([2,3]) );
      C.push_back( XY([3,4]) );
      C.push_back( XY([4,5]) );      
      C.push_back( XY([5,6]) );
      
      % layer 4
      C.push_back( XY([1,2]) );
      C.push_back( XY([2,3]) );
      C.push_back( XY([3,4]) );
      C.push_back( XY([4,5]) );      
      C.push_back( XY([5,6]) );
      
      % layer 5            
      C.push_back( XY([0,1]) );
      C.push_back( XY([1,2]) );
      C.push_back( XY([2,3]) );
      C.push_back( XY([3,4]) );
      C.push_back( XY([4,5]) );      
      C.push_back( XY([5,6]) );
      
      ascend = false;
      DTrC = TrCircuit( C, ascend );
      
      % check index functions
      for i = 1:DTrC.nbGates
        [q, l] = DTrC.QLIndex( i );
        test.verifyEqual( DTrC.descendLinearIndex( q, l ), i );
      end
      
      test.verifyEqual( DTrC.ascendLinearIndex( 4, 0 ), 2 );
      
      % ascend / descend
      test.verifyFalse( DTrC.ascend );
      test.verifyTrue( DTrC.descend );
    end
    
    function test_TriangleCircuit_Ascend_Descend_reorder( test )
      tol = 100 * eps;
      TrC = test.random7QubitAscendTriangle;
      
      test.verifyTrue( TrC.ascend );
      test.verifyFalse( TrC.descend );
      
      matAscend = TrC.matrix;
      
      TrC.makeDescend();
      
      test.verifyFalse( TrC.ascend );
      test.verifyTrue( TrC.descend );
      
      test.verifyEqual( TrC.matrix, matAscend, 'AbsTol', tol );
      
      TrC.makeAscend();
      
      test.verifyTrue( TrC.ascend );
      test.verifyFalse( TrC.descend );
      
      test.verifyEqual( TrC.matrix, matAscend ); % should be exactly equal
      
      TrC = test.random8QubitDescendTriangle;
      
      test.verifyFalse( TrC.ascend );
      test.verifyTrue( TrC.descend );
      
      matDescend = TrC.matrix;
      
      TrC.makeAscend();
      
      test.verifyTrue( TrC.ascend );
      test.verifyFalse( TrC.descend );
      
      test.verifyEqual( TrC.matrix, matDescend, 'AbsTol', tol );
      
      TrC.makeDescend();
      
      test.verifyFalse( TrC.ascend );
      test.verifyTrue( TrC.descend );
      
      test.verifyEqual( TrC.matrix, matDescend ); % should be exactly equal
    end
    
    function test_Replace( test )
      TrC = test.random7QubitAscendTriangle;
      
      XY = @f3c.qgates.RotationXY ;
      gate = XY( [5,6], 0, 0 );
      
      gate1 = TrC.gates(1);
      test.verifyNotEqual( gate1.thetas, [0, 0]);
      
      TrC.replace( 1, gate );
      gate1 = TrC.gates(1);
      test.verifyEqual( gate1.qubits, gate.qubits );
      test.verifyEqual( gate1.thetas, [0, 0] );
      
    end
    
    function test_TriangleCircuit_mergeGate( test )
      tol = 100 * eps;
      XY = @f3c.qgates.RotationXY ;
      Circuit = @qclab.QCircuit ;
      
      % Ascending 7-qubit Triangle  
      C = Circuit( 7 );
      
      % left merge -- lowest qubit, no turnovers
      ATrC = test.random7QubitAscendTriangle;      
      ATrCmat = ATrC.matrix ;
      
      gate = XY([5,6], randn, randn );
      C.push_back( gate );
      Cmat = C.matrix ;
      
      ATrC.merge( gate, 'L' );
      test.verifyEqual( ATrC.matrix, ATrCmat * Cmat, 'AbsTol', tol );
      
      % right merge -- lowest qubit, no turnovers
      ATrC = test.random7QubitAscendTriangle;      
      ATrCmat = ATrC.matrix ;
      
      gate = XY([5,6], randn, randn );
      C.clear();
      C.push_back( gate );
      Cmat = C.matrix ;
      
      ATrC.merge( gate, 'R' );
      test.verifyEqual( ATrC.matrix, Cmat * ATrCmat, 'AbsTol', tol );
      
      % left merge -- first qubit, maximum turnovers
      ATrC = test.random7QubitAscendTriangle;      
      ATrCmat = ATrC.matrix ;
      
      gate = XY([0,1], randn, randn );
      C.clear();
      C.push_back( gate );
      Cmat = C.matrix ;
      
      ATrC.merge( gate, 'L' );
      test.verifyEqual( ATrC.matrix, ATrCmat * Cmat, 'AbsTol', tol );
      
      % right merge -- first qubit, maximum turnovers
      ATrC = test.random7QubitAscendTriangle;      
      ATrCmat = ATrC.matrix ;
      
      gate = XY([0,1], randn, randn );
      C.clear();
      C.push_back( gate );
      Cmat = C.matrix ;
      
      ATrC.merge( gate, 'R' );
      test.verifyEqual( ATrC.matrix, Cmat * ATrCmat, 'AbsTol', tol );
      
      % left merge -- qubit in the middle
       ATrC = test.random7QubitAscendTriangle;      
      ATrCmat = ATrC.matrix ;
      
      gate = XY([2,3], randn, randn );
      C.clear();
      C.push_back( gate );
      Cmat = C.matrix ;
      
      ATrC.merge( gate, 'L' );
      test.verifyEqual( ATrC.matrix, ATrCmat * Cmat, 'AbsTol', tol );
      
      % right merge -- qubit in the middle
      ATrC = test.random7QubitAscendTriangle;      
      ATrCmat = ATrC.matrix ;
      
      gate = XY([3,4], randn, randn );
      C.clear();
      C.push_back( gate );
      Cmat = C.matrix ;
      
      ATrC.merge( gate, 'R' );
      test.verifyEqual( ATrC.matrix, Cmat * ATrCmat, 'AbsTol', tol );
      
      % Descending 8-qubit Triangle  
      C = Circuit( 8 );
      
      % left merge -- lowest qubit, no turnovers
      DTrC = test.random8QubitDescendTriangle;      
      DTrCmat = DTrC.matrix ;
      
      gate = XY([6, 7], randn, randn );
      C.clear();
      C.push_back( gate );
      Cmat = C.matrix ;
      
      DTrC.merge( gate, 'L' );
      test.verifyEqual( DTrC.matrix, DTrCmat * Cmat, 'AbsTol', tol );
      
      % right merge - lowest qubit, no turnovers
      DTrC = test.random8QubitDescendTriangle;      
      DTrCmat = DTrC.matrix ;
      
      gate = XY([6, 7], randn, randn );
      C.clear();
      C.push_back( gate );
      Cmat = C.matrix ;
      
      DTrC.merge( gate, 'R' );
      test.verifyEqual( DTrC.matrix, Cmat * DTrCmat, 'AbsTol', tol );
      
      % left merge -- first qubit, maximum turnovers
      DTrC = test.random8QubitDescendTriangle;      
      DTrCmat = DTrC.matrix ;
      
      gate = XY([0, 1], randn, randn );
      C.clear();
      C.push_back( gate );
      Cmat = C.matrix ;
      
      DTrC.merge( gate, 'L' );
      test.verifyEqual( DTrC.matrix, DTrCmat * Cmat, 'AbsTol', tol );
      
      % right merge - first qubit, maximum turnovers
      DTrC = test.random8QubitDescendTriangle;      
      DTrCmat = DTrC.matrix ;
      
      gate = XY([0, 1], randn, randn );
      C.clear();
      C.push_back( gate );
      Cmat = C.matrix ;
      
      DTrC.merge( gate, 'R' );
      test.verifyEqual( DTrC.matrix, Cmat * DTrCmat, 'AbsTol', tol );
      
      % left merge -- qubit in the middle
      DTrC = test.random8QubitDescendTriangle;      
      DTrCmat = DTrC.matrix ;
      
      gate = XY([3,4], randn, randn );
      C.clear();
      C.push_back( gate );
      Cmat = C.matrix ;
      
      DTrC.merge( gate, 'L' );
      test.verifyEqual( DTrC.matrix, DTrCmat * Cmat, 'AbsTol', tol );
      
      % right merge - qubit in the middle
      DTrC = test.random8QubitDescendTriangle;      
      DTrCmat = DTrC.matrix ;
      
      gate = XY([4, 5], randn, randn );
      C.clear();
      C.push_back( gate );
      Cmat = C.matrix ;
      
      DTrC.merge( gate, 'R' );
      test.verifyEqual( DTrC.matrix, Cmat * DTrCmat, 'AbsTol', tol );
      
    end
    
    function test_TriangleCircuit_to_SquareCircuit( test )
      tol = 1000 * eps;
      
      % Ascending Triangle on odd number of qubits
      ATrC = test.random7QubitAscendTriangle;
      ATrCmat = ATrC.matrix ;
      
      % Transform to Square
      SqC = f3c.SquareCircuit( ATrC );
      
      test.verifyEqual( SqC.matrix, ATrCmat, 'AbsTol', tol );
      test.verifyEqual( ATrC.matrix, ATrCmat ); % should be exactly equal
      
      % Descending Triangle on even number of qubits
      DTrC = test.random8QubitDescendTriangle;
      DTrCmat = DTrC.matrix;
      
      % Transform to Square
      SqC = f3c.SquareCircuit( DTrC );
      
      test.verifyEqual( SqC.matrix, DTrCmat, 'AbsTol', tol );
      test.verifyEqual( DTrC.matrix, DTrCmat ); % should be exactly equal
    end
  end
  
  methods (Static)
    function ATrC = random7QubitAscendTriangle()
      Circuit = @qclab.QCircuit ;
      TrCircuit = @f3c.TriangleCircuit ;
      XY = @f3c.qgates.RotationXY ;
      
      % Ascending Triangle
      C = Circuit( 7 );
      
      % layer 0
      C.push_back( XY([5,6], randn, randn ) );
      C.push_back( XY([4,5], randn, randn ) );      
      C.push_back( XY([3,4], randn, randn ) );
      C.push_back( XY([2,3], randn, randn ) );
      C.push_back( XY([1,2], randn, randn ) );
      C.push_back( XY([0,1], randn, randn ) );
      
      % layer 1      
      C.push_back( XY([5,6], randn, randn ) );
      C.push_back( XY([4,5], randn, randn ) );      
      C.push_back( XY([3,4], randn, randn ) );
      C.push_back( XY([2,3], randn, randn ) );
      C.push_back( XY([1,2], randn, randn ) );
      
      % layer 2
      C.push_back( XY([5,6], randn, randn ) );
      C.push_back( XY([4,5], randn, randn ) );      
      C.push_back( XY([3,4], randn, randn ) );
      C.push_back( XY([2,3], randn, randn ) );
      
      % layer 3
      C.push_back( XY([5,6], randn, randn ) );
      C.push_back( XY([4,5], randn, randn ) );      
      C.push_back( XY([3,4], randn, randn ) );
      
      % layer 4
      C.push_back( XY([5,6], randn, randn ) );
      C.push_back( XY([4,5], randn, randn ) );      
      
      % layer 5
      C.push_back( XY([5,6], randn, randn ) );
      
      ascend = true;
      ATrC = TrCircuit( C, ascend );
    end
    
    function DTrC = random8QubitDescendTriangle()
      Circuit = @qclab.QCircuit ;
      TrCircuit = @f3c.TriangleCircuit ;
      XY = @f3c.qgates.RotationXY ;
      
      % Ascending Triangle
      C = Circuit( 8 );
      
      % layer 0
      C.push_back( XY([6,7], randn, randn ) );
      
      % layer 1
      C.push_back( XY([5,6], randn, randn ) );
      C.push_back( XY([6,7], randn, randn ) );
      
      % layer 2
      C.push_back( XY([4,5], randn, randn ) );      
      C.push_back( XY([5,6], randn, randn ) );
      C.push_back( XY([6,7], randn, randn ) );
      
      % layer 3
      C.push_back( XY([3,4], randn, randn ) );
      C.push_back( XY([4,5], randn, randn ) );      
      C.push_back( XY([5,6], randn, randn ) );
      C.push_back( XY([6,7], randn, randn ) );
      
      % layer 4
      C.push_back( XY([2,3], randn, randn ) );
      C.push_back( XY([3,4], randn, randn ) );
      C.push_back( XY([4,5], randn, randn ) );      
      C.push_back( XY([5,6], randn, randn ) );
      C.push_back( XY([6,7], randn, randn ) );
      
      % layer 5
      C.push_back( XY([1,2], randn, randn ) );
      C.push_back( XY([2,3], randn, randn ) );
      C.push_back( XY([3,4], randn, randn ) );
      C.push_back( XY([4,5], randn, randn ) );      
      C.push_back( XY([5,6], randn, randn ) );
      C.push_back( XY([6,7], randn, randn ) );
      
      % layer 6            
      C.push_back( XY([0,1], randn, randn ) );
      C.push_back( XY([1,2], randn, randn ) );
      C.push_back( XY([2,3], randn, randn ) );
      C.push_back( XY([3,4], randn, randn ) );
      C.push_back( XY([4,5], randn, randn ) );      
      C.push_back( XY([5,6], randn, randn ) );
      C.push_back( XY([6,7], randn, randn ) );
      
      ascend = false;
      DTrC = TrCircuit( C, ascend );
    end
    
  end
end