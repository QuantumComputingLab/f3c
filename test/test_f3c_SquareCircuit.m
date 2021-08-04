classdef test_f3c_SquareCircuit < matlab.unittest.TestCase
  methods (Test)
    function test_SquareCircuit( test)
      Circuit = @qclab.QCircuit ;
      SqCircuit = @f3c.SquareCircuit ;
      XY = @f3c.qgates.RotationXY ;
      
      % Even number of qubits
      % Generate circuit
      C = Circuit( 4 );
      C.push_back( XY([0,1]) );
      C.push_back( XY([2,3]) );
      C.push_back( XY([1,2]) );      
      C.push_back( XY([0,1]) );
      C.push_back( XY([2,3]) );
      C.push_back( XY([1,2]) );
      
      SqC = SqCircuit( C );
      
      % qubits
      test.verifyEqual( SqC.nbQubits, int32(4) );
      test.verifyEqual( SqC.qubit, int32(0) );
      test.verifyEqual( SqC.qubits, int32(0:3) );
      
      % nbGates
      test.verifyEqual( SqC.nbGates, 6 );
      
      % offset
      test.verifyEqual( SqC.offset, int32(0) );
      
      % matrix
      test.verifyEqual( SqC.matrix, qclab.qId(4), 'AbsTol', eps(100) );
      
      % replace test
      
      % Check that SquareCircuit can't be modified
      % setOffset
      myHandle = @(in, in2) f3c.SquareCircuit.setOffset(in, in2) ;
      test.verifyError( @() myHandle(SqC, XY([7,10])), 'MATLAB:assertion:failed');
      % gateHandles
      myHandle = @(in, in2) f3c.SquareCircuit.gateHandles(in, in2) ;
      test.verifyError( @() myHandle(SqC, XY([7,10])), 'MATLAB:assertion:failed');
      % clear
      myHandle = @(in) f3c.SquareCircuit.clear(in) ;
      test.verifyError( @() myHandle(SqC), 'MATLAB:assertion:failed' );
      % erase
      myHandle = @(in) f3c.SquareCircuit.erase(in) ;
      test.verifyError( @() myHandle(SqC), 'MATLAB:assertion:failed' );
      % push_back
      myHandle = @(in, in2) f3c.SquareCircuit.push_back(in, in2) ;
      test.verifyError( @() myHandle(SqC, XY([7,10])), 'MATLAB:assertion:failed');
      % pop_back
      myHandle = @(in) f3c.SquareCircuit.pop_back(in) ;
      test.verifyError( @() myHandle(SqC), 'MATLAB:assertion:failed' );
      
      % Odd number of qubits
      C = Circuit( 5 );
      % layer 0
      C.push_back( XY([0,1]) );
      C.push_back( XY([2,3]) );
      
      % layer 1
      C.push_back( XY([1,2]) );      
      C.push_back( XY([3,4]) );
      
      % layer 2
      C.push_back( XY([0,1]) );
      C.push_back( XY([2,3]) );
      
      % layer 3
      C.push_back( XY([1,2]) );
      C.push_back( XY([3,4]) );
      
      % layer 4
      C.push_back( XY([0,1]) );
      C.push_back( XY([2,3]) );
      
      SqC = SqCircuit( C );
      
      % qubits
      test.verifyEqual( SqC.nbQubits, int32(5) );
      test.verifyEqual( SqC.qubit, int32(0) );
      test.verifyEqual( SqC.qubits, int32(0:4) );
      
      % nbGates
      test.verifyEqual( SqC.nbGates, 10 );
      
      % check index functions
      for i = 1:SqC.nbGates
        [q, l] = SqC.qlIndex( i );
        test.verifyEqual( SqC.linearIndex( q, l ), i );
      end
      
      C = Circuit( 7 );
      
      % layer 0
      C.push_back( XY([0,1]) );
      C.push_back( XY([2,3]) );
      C.push_back( XY([4,5]) );
      
      % layer 1
      C.push_back( XY([1,2]) );      
      C.push_back( XY([3,4]) );
      C.push_back( XY([5,6]) );
      
      % layer 2
      C.push_back( XY([0,1]) );
      C.push_back( XY([2,3]) );
      C.push_back( XY([4,5]) );
      
      % layer 3
      C.push_back( XY([1,2]) );      
      C.push_back( XY([3,4]) );
      C.push_back( XY([5,6]) );
      
      % layer 4
      C.push_back( XY([0,1]) );
      C.push_back( XY([2,3]) );
      C.push_back( XY([4,5]) );
      
      % layer 5
      C.push_back( XY([1,2]) );      
      C.push_back( XY([3,4]) );
      C.push_back( XY([5,6]) );
      
      % layer 6
      C.push_back( XY([0,1]) );
      C.push_back( XY([2,3]) );
      C.push_back( XY([4,5]) );
      
      
      SqC = SqCircuit( C );
      
      % check index functions
      for i = 1:SqC.nbGates
        [q, l] = SqC.qlIndex( i );
        test.verifyEqual( SqC.linearIndex( q, l ), i );
      end
      
      C = Circuit( 8 );
      
      % layer 0
      C.push_back( XY([0,1]) );
      C.push_back( XY([2,3]) );
      C.push_back( XY([4,5]) );
      C.push_back( XY([6,7]) );
      
      % layer 1
      C.push_back( XY([1,2]) );      
      C.push_back( XY([3,4]) );
      C.push_back( XY([5,6]) );
      
      % layer 2
      C.push_back( XY([0,1]) );
      C.push_back( XY([2,3]) );
      C.push_back( XY([4,5]) );
      C.push_back( XY([6,7]) );
      
      % layer 3
      C.push_back( XY([1,2]) );      
      C.push_back( XY([3,4]) );
      C.push_back( XY([5,6]) );
      
      % layer 4
      C.push_back( XY([0,1]) );
      C.push_back( XY([2,3]) );
      C.push_back( XY([4,5]) );
      C.push_back( XY([6,7]) );
      
      % layer 5
      C.push_back( XY([1,2]) );      
      C.push_back( XY([3,4]) );
      C.push_back( XY([5,6]) );
      
      % layer 6
      C.push_back( XY([0,1]) );
      C.push_back( XY([2,3]) );
      C.push_back( XY([4,5]) );
      C.push_back( XY([6,7]) );
      
      % layer 7
      C.push_back( XY([1,2]) );      
      C.push_back( XY([3,4]) );
      C.push_back( XY([5,6]) );
      
      
      SqC = SqCircuit( C );
      
      % check index functions
      for i = 1:SqC.nbGates
        [q, l] = SqC.qlIndex( i );
        test.verifyEqual( SqC.linearIndex( q, l ), i );
      end
      
    end
    
    function test_Replace( test )
      SqC = test.random8QubitSquare;
      
      XY = @f3c.qgates.RotationXY ;
      gate = XY( [0,1], 0, 0 );
      
      gate1 = SqC.gates(1);
      test.verifyNotEqual( gate1.thetas, [0, 0]);
      
      SqC.replace( 1, gate );
      gate1 = SqC.gates(1);
      test.verifyEqual( gate1.qubits, gate.qubits );
      test.verifyEqual( gate1.thetas, [0, 0] );
      
    end
    
    function test_Square_to_Triangle( test )
      tol = 1000 * eps;
      
      % Square on even qubits
      SqC = test.random8QubitSquare;
      SqCmat = SqC.matrix ;
      
      % Transform to Triangle
      TrC = f3c.TriangleCircuit( SqC );
      
      test.verifyEqual( TrC.matrix, SqCmat, 'AbsTol', tol );
      test.verifyEqual( SqC.matrix, SqCmat ); % should be exactly equal
      
      % Square on odd qubits
      SqC = test.random9QubitSquare;
      SqCmat = SqC.matrix ;
      
      % Transform to Triangle
      TrC = f3c.TriangleCircuit( SqC );
      
      test.verifyEqual( TrC.matrix, SqCmat, 'AbsTol', tol );
      test.verifyEqual( SqC.matrix, SqCmat ); % should be exactly equal
    end
  end
  
  methods (Static)
    function SqC = random8QubitSquare()
      Circuit = @qclab.QCircuit ;
      SqCircuit = @f3c.SquareCircuit ;
      XY = @f3c.qgates.RotationXY ;
      
      % Square
      C = Circuit( 8 );
      
      % layer 0            
      C.push_back( XY([0,1], randn, randn ) );
      C.push_back( XY([2,3], randn, randn ) );
      C.push_back( XY([4,5], randn, randn ) );
      C.push_back( XY([6,7], randn, randn ) );
      
      % layer 1
      C.push_back( XY([1,2], randn, randn ) );
      C.push_back( XY([3,4], randn, randn ) );
      C.push_back( XY([5,6], randn, randn ) );
      
      % layer 2
      C.push_back( XY([0,1], randn, randn ) );
      C.push_back( XY([2,3], randn, randn ) );
      C.push_back( XY([4,5], randn, randn ) );
      C.push_back( XY([6,7], randn, randn ) );
      
      % layer 3
      C.push_back( XY([1,2], randn, randn ) );
      C.push_back( XY([3,4], randn, randn ) );
      C.push_back( XY([5,6], randn, randn ) );
      
      % layer 4
      C.push_back( XY([0,1], randn, randn ) );
      C.push_back( XY([2,3], randn, randn ) );
      C.push_back( XY([4,5], randn, randn ) );
      C.push_back( XY([6,7], randn, randn ) );
      
      % layer 5
      C.push_back( XY([1,2], randn, randn ) );
      C.push_back( XY([3,4], randn, randn ) );
      C.push_back( XY([5,6], randn, randn ) );
      
      % layer 6
      C.push_back( XY([0,1], randn, randn ) );
      C.push_back( XY([2,3], randn, randn ) );
      C.push_back( XY([4,5], randn, randn ) );
      C.push_back( XY([6,7], randn, randn ) );
      
      % layer 7
      C.push_back( XY([1,2], randn, randn ) );
      C.push_back( XY([3,4], randn, randn ) );
      C.push_back( XY([5,6], randn, randn ) );
      
      SqC = SqCircuit( C );
    end
    
    function SqC = random9QubitSquare()
      Circuit = @qclab.QCircuit ;
      SqCircuit = @f3c.SquareCircuit ;
      XY = @f3c.qgates.RotationXY ;
      
      % Square
      C = Circuit( 9 );
      
      % layer 0            
      C.push_back( XY([0,1], randn, randn ) );
      C.push_back( XY([2,3], randn, randn ) );
      C.push_back( XY([4,5], randn, randn ) );
      C.push_back( XY([6,7], randn, randn ) );
      
      % layer 1
      C.push_back( XY([1,2], randn, randn ) );
      C.push_back( XY([3,4], randn, randn ) );
      C.push_back( XY([5,6], randn, randn ) );
      C.push_back( XY([7,8], randn, randn ) );     
      
      % layer 2            
      C.push_back( XY([0,1], randn, randn ) );
      C.push_back( XY([2,3], randn, randn ) );
      C.push_back( XY([4,5], randn, randn ) );
      C.push_back( XY([6,7], randn, randn ) );
      
      % layer 3
      C.push_back( XY([1,2], randn, randn ) );
      C.push_back( XY([3,4], randn, randn ) );
      C.push_back( XY([5,6], randn, randn ) );
      C.push_back( XY([7,8], randn, randn ) );
         
      % layer 4            
      C.push_back( XY([0,1], randn, randn ) );
      C.push_back( XY([2,3], randn, randn ) );
      C.push_back( XY([4,5], randn, randn ) );
      C.push_back( XY([6,7], randn, randn ) );
      
      % layer 5
      C.push_back( XY([1,2], randn, randn ) );
      C.push_back( XY([3,4], randn, randn ) );
      C.push_back( XY([5,6], randn, randn ) );
      C.push_back( XY([7,8], randn, randn ) );
        
      % layer 6            
      C.push_back( XY([0,1], randn, randn ) );
      C.push_back( XY([2,3], randn, randn ) );
      C.push_back( XY([4,5], randn, randn ) );
      C.push_back( XY([6,7], randn, randn ) );
      
      % layer 7
      C.push_back( XY([1,2], randn, randn ) );
      C.push_back( XY([3,4], randn, randn ) );
      C.push_back( XY([5,6], randn, randn ) );
      C.push_back( XY([7,8], randn, randn ) );
      
      % layer 8            
      C.push_back( XY([0,1], randn, randn ) );
      C.push_back( XY([2,3], randn, randn ) );
      C.push_back( XY([4,5], randn, randn ) );
      C.push_back( XY([6,7], randn, randn ) );
      
      SqC = SqCircuit( C );
    end
  end
end