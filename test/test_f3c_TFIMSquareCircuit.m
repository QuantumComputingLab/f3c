classdef test_f3c_TFIMSquareCircuit < matlab.unittest.TestCase
  methods (Test)
    function test_TFIMSquareCircuit( test)
      
      % Create a TFIM Square Circuit on 4 spins.
      
      C = test.randomTFIMCircuit(4);      
      SqC = f3c.TFIMSquareCircuit( C );
      
      % qubits
      test.verifyEqual( SqC.nbQubits, int32(4) );
      test.verifyEqual( SqC.qubit, int32(0) );
      test.verifyEqual( SqC.qubits, int32(0:3) );
      
      % nbGates
      test.verifyEqual( SqC.nbGates, 28 );
      
      % offset
      test.verifyEqual( SqC.offset, int32(0) );
      
      % matrix
      SqCmat = SqC.matrix;
      test.verifyEqual( SqCmat, C.matrix, 'AbsTol', eps(100) );
      
      % check index functions
      for i = 1:SqC.nbGates
        [l, j] = SqC.linearTo2DIndex( SqC.nbQubits, i );
        test.verifyEqual( SqC.linearIndex( SqC.nbQubits, l, j ), i );
      end
      
      % transform to an ascending Square circuit
      ATrC = f3c.TFIMTriangleCircuit( SqC );
      test.verifyEqual( SqCmat, ATrC.matrix, 'AbsTol', eps(100) );
      
    end
  end
  
  methods (Static)
    function [C] = randomTFIMCircuit(n)
      Z = @qclab.qgates.RotationZ;
      XX = @qclab.qgates.RotationXX;
      
      C = qclab.QCircuit( n );
      
      for l = 0:n-1
        % single spins
        for j = 0:n-1
          C.push_back( Z(j, randn) );
        end
        % two spins
        for j=0:n-2
          C.push_back( XX([j,j+1], randn) );
        end
      end
      
    end
  end
end