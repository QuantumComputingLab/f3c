classdef test_f3c_IsingCircuit < matlab.unittest.TestCase
  methods (Test)
    function test_IsingCircuit_basic( test )
      IsingCirc = f3c.IsingCircuit( 2, 'X' );
      FullCirc = qclab.QCircuit( 2 );

      IsingCirc.mergeOneQRot( qclab.qgates.RotationX( 0, 0.1 ) );
      FullCirc.push_back( qclab.qgates.RotationX( 0, 0.1 ) );
      test.verifyEqual( IsingCirc.matrix, FullCirc.matrix, 'AbsTol', 10*eps ) ;
      
      IsingCirc.mergeOneQRot( qclab.qgates.RotationX( 1, 0.1 ) );
      FullCirc.push_back( qclab.qgates.RotationX( 1, 0.1 ) );
      test.verifyEqual( IsingCirc.matrix, FullCirc.matrix, 'AbsTol', 10*eps ) ;
      
      IsingCirc.mergeOneQRot( qclab.qgates.RotationX( 0, 0.2 ) );
      FullCirc.push_back( qclab.qgates.RotationX( 0, 0.2 ) );
      test.verifyEqual( IsingCirc.matrix, FullCirc.matrix, 'AbsTol', 10*eps ) ;
      
      IsingCirc.mergeTwoQRot( qclab.qgates.RotationXX( [0, 1], 0.3 ) );
      FullCirc.push_back( qclab.qgates.RotationXX( [0, 1], 0.3 ) );
      test.verifyEqual( IsingCirc.matrix, FullCirc.matrix, 'AbsTol', 10*eps ) ;
      
      IsingCirc.mergeOneQRot( qclab.qgates.RotationX( 1, 0.4 ) );
      FullCirc.push_back( qclab.qgates.RotationX( 1, 0.4 ) );
      test.verifyEqual( IsingCirc.matrix, FullCirc.matrix, 'AbsTol', 10*eps ) ;
      
      IsingCirc.mergeTwoQRot( qclab.qgates.RotationXX( [0, 1], 0.5 ) );
      FullCirc.push_back( qclab.qgates.RotationXX( [0, 1], 0.5 ) );
      test.verifyEqual( IsingCirc.matrix, FullCirc.matrix, 'AbsTol', 10*eps ) ;
    end
  end
end

