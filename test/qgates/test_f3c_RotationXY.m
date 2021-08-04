classdef test_f3c_RotationXY < matlab.unittest.TestCase
  methods (Test)
    function test_RotationXY_basic_functionality(test)
      XY = @f3c.qgates.RotationXY ;
      
      % default constructor
      G1 = XY() ;
      test.verifyEqual( G1.nbQubits, int32(2) );    % nbQubits
      test.verifyFalse( G1.fixed );                 % fixed
      test.verifyFalse( G1.controlled );            % controlled
      
      % matrix
      test.verifyEqual( G1.matrix, eye(4) );
      
      % qubits
      qubits = G1.qubits;
      test.verifyEqual( length(qubits), 2 );
      test.verifyEqual( qubits(1), int32(0) );
      test.verifyEqual( qubits(2), int32(1) );
      qnew = [5, 3] ;
      G1.setQubits( qnew );
      qubits = G1.qubits;
      test.verifyEqual( qubits(1), int32(3) );
      test.verifyEqual( qubits(2), int32(5) );
      
      % fixed
      G1.makeFixed();
      test.verifyTrue( G1.fixed );
      G1.makeVariable();
      test.verifyFalse( G1.fixed );
    end
    
    function test_RotationXY_constructors(test)
      tol = 10*eps;
      XY = @f3c.qgates.RotationXY ;
      
      % 2.
      G2 = XY([4,2]);
      qubits = G2.qubits;
      test.verifyEqual( length(qubits), 2 );
      test.verifyEqual( qubits(1), int32(2) );
      test.verifyEqual( qubits(2), int32(4) );
      
      % 3.
      theta1 = -1.2;
      theta2 = 2.2;
      
      G3 = XY( [4,2], theta1, theta2 );
      
      qubits = G3.qubits;
      test.verifyEqual( length(qubits), 2 );
      test.verifyEqual( qubits(1), int32(2) );
      test.verifyEqual( qubits(2), int32(4) );
      
      % circuit angles
      thetas = G3.thetas ;
      modDif = mod(thetas(1) - theta1, 2*pi);
      test.verifyTrue( min(modDif, modDif - 2*pi) < tol  );
      modDif = mod(thetas(2) - theta2, 2*pi);
      test.verifyTrue( min(modDif, modDif - 2*pi) < tol  );
      
      % 4.
      rot1 = qclab.QRotation( theta1 );
      rot2 = qclab.QRotation( theta2 );
      G4 = XY( [4,2], rot1, rot2 );
      test.verifyTrue( G3 == G4 );
    end
    
    function test_RotationXY_update(test)
      tol = 10*eps;
      XY = @f3c.qgates.RotationXY ;
      
      G1 = XY() ;
      
      % 1.
      theta1 = pi/4;
      theta2 = -pi/4;
      G1.update( theta1, theta2 );
      
      % circuit angles
      thetas = G1.thetas ;
      modDif = mod(thetas(1) - theta1, 2*pi);
      test.verifyTrue( min(modDif, modDif - 2*pi) < tol  );
      modDif = mod(thetas(2) - theta2, 2*pi);
      test.verifyTrue( min(modDif, modDif - 2*pi) < tol  );
      
      % 2.
      rot1 = qclab.QRotation( theta1 );
      rot2 = qclab.QRotation( theta2 );
      
      G2 = XY() ;
      G2.update( rot1, rot2 );
      test.verifyTrue( G1 == G2 );
    end
    
    function test_RotationXY_matrix( test )
      tol = 10*eps;
      XY = @f3c.qgates.RotationXY ;
      
      G1 = XY() ;
      theta1 = pi/7;
      theta2 = 5.33;
      G1.update( theta1, theta2 );
      cplus = cos( (theta1 + theta2)/2);
      splus = sin( (theta1 + theta2)/2);
      cminus = cos( (theta1 - theta2)/2);
      sminus = sin( (theta1 - theta2)/2);
      
      check = [cminus, 0, 0, -1i*sminus;
               0, cplus, -1i*splus, 0;
               0, -1i*splus, cplus, 0;
               -1i*sminus, 0, 0, cminus];
      test.verifyEqual( G1.matrix, check, 'AbsTol', tol );
    end
    
    function test_RotationXY_fuse_multiply( test )
      tol = 10*eps;
      XY = @f3c.qgates.RotationXY ;
      
      G1 = XY() ;
      G1.update( pi/4, -pi/4 );
      G2 = XY() ;
      G2.update( 8*pi/4, pi/7 );
      G1c = copy(G1);
      G1.fuse( G2, 'R' );
      test.verifyEqual( G1.matrix, G2.matrix * G1c.matrix, 'AbsTol', tol );
      G1 = copy(G1c);
      G1.fuse( G2, 'L' );
      test.verifyEqual( G1.matrix, G1c.matrix * G2.matrix, 'AbsTol', tol );
      
      
      G12 = G1c * G2;
      test.verifyEqual( G1.matrix, G12.matrix , 'AbsTol', tol );
    end
    
    function  test_RotationXY_turnover( test )
      tol = 10*eps;
      XY = @f3c.qgates.RotationXY ;
      Circuit = @qclab.QCircuit ;
      
      % test Vee to Hat direction
      theta1 = [ pi/7, 5.33];
      G1 = XY([0, 1], theta1(1), theta1(2) );
      
      theta2 = [ pi/4, -2.21];
      G2 = XY([1, 2], theta2(1), theta2(2) );
      
      theta3 = [ -0.81, pi/5];
      G3 = XY([0, 1], theta3(1), theta3(2) );
      
      C = Circuit(3);
      C.push_back( G1 );
      C.push_back( G2 );
      C.push_back( G3 );
      
      Mat = C.matrix ;
      
      [GA, GB, GC] = turnover( G1, G2, G3 );

      C.replace( 1, GA );
      C.replace( 2, GB );
      C.replace( 3, GC );
      
      % check qubits after turnover
      qubits1 = G1.qubits;
      qubits2 = G2.qubits;
      qubitsA = GA.qubits;
      test.verifyEqual( qubitsA, qubits2 );
      qubitsB = GB.qubits;
      test.verifyEqual( qubitsB, qubits1 );
      qubitsC = GC.qubits;
      test.verifyEqual( qubitsC, qubits2 );
      
      % check circuit after turnover and compare with before
      test.verifyEqual( C.matrix, Mat, 'AbsTol', tol ) ;
      
      % test Hat to Vee direction
      G1 = XY([1, 2], theta1(1), theta1(2) );      
      G2 = XY([0, 1], theta2(1), theta2(2) );      
      G3 = XY([1, 2], theta3(1), theta3(2) );
      
      C = Circuit(3);
      C.push_back( G1 );
      C.push_back( G2 );
      C.push_back( G3 );
      
      Mat = C.matrix ;
      
      fprintf(1, '\n');
      C.draw( 1, 'S' );
      [GA, GB, GC] = turnover( G1, G2, G3 );
      C.replace( 1, GA );
      C.replace( 2, GB );
      C.replace( 3, GC );
      C.draw( 1, 'S' );
      fprintf(1, '\n');

      % check qubits after turnover
      qubits1 = G1.qubits;
      qubits2 = G2.qubits;
      qubitsA = GA.qubits;
      test.verifyEqual( qubitsA, qubits2 );
      qubitsB = GB.qubits;
      test.verifyEqual( qubitsB, qubits1 );
      qubitsC = GC.qubits;
      test.verifyEqual( qubitsC, qubits2 );
      
      % check circuit after turnover and compare with before
      test.verifyEqual( C.matrix, Mat, 'AbsTol', tol ) ;
    end
  end
end