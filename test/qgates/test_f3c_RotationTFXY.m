classdef test_f3c_RotationTFXY < matlab.unittest.TestCase
  methods (Test)
    function test_RotationTFXY_basic_functionality(test)
      TFXY = @f3c.qgates.RotationTFXY ;
      
      % default constructor
      G1 = TFXY() ;
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
    
    function test_RotationTFXY_constructors(test)
      tol = 10*eps;
      TFXY = @f3c.qgates.RotationTFXY ;
      
      % 2.
      G2 = TFXY([4,2]);
      qubits = G2.qubits;
      test.verifyEqual( length(qubits), 2 );
      test.verifyEqual( qubits(1), int32(2) );
      test.verifyEqual( qubits(2), int32(4) );
      
      % 3.
      theta1 = -1.2;
      theta2 = 2.2;
      theta3 = 0.75;
      theta4 = -2.11;
      theta5 = -3.09;
      theta6 = -0.88;
      
      G3 = TFXY( [4,2], theta1, theta2, theta3, theta4, theta5, theta6 );
      
      qubits = G3.qubits;
      test.verifyEqual( length(qubits), 2 );
      test.verifyEqual( qubits(1), int32(2) );
      test.verifyEqual( qubits(2), int32(4) );
      
      % left circuit angles
      lTheta = G3.thetasLeft ;
      modDif = mod(lTheta(1) - theta1 - theta2, 2*pi);
      test.verifyTrue( min(modDif, modDif - 2*pi) < tol  );
      modDif = mod(lTheta(2) - theta1 + theta2, 2*pi);
      test.verifyTrue( min(modDif, modDif - 2*pi) < tol  );
      % middle circuit angles
      mTheta = G3.thetasMiddle ;
      modDif = mod(mTheta(1) - theta3 - theta4, 2*pi);
      test.verifyTrue( min(modDif, modDif - 2*pi) < tol  );
      modDif = mod(mTheta(2) - theta3 + theta4, 2*pi);
      test.verifyTrue( min(modDif, modDif - 2*pi) < tol  );
      % right circuit angles
      rTheta = G3.thetasRight ;
      modDif = mod(rTheta(1) - theta5 - theta6, 2*pi);
      test.verifyTrue( min(modDif, modDif - 2*pi) < tol  );
      modDif = mod(rTheta(2) - theta5 + theta6, 2*pi);
      test.verifyTrue( min(modDif, modDif - 2*pi) < tol  );
      
      % 4.
      theta = [ theta1, theta2, theta3, theta4, theta5, theta6 ];
      G4 = TFXY( [4,2], theta );
      test.verifyTrue( G3 == G4 );
      
      % 5.
      rot1 = qclab.QRotation( theta1 );
      rot2 = qclab.QRotation( theta2 );
      rot3 = qclab.QRotation( theta3 );
      rot4 = qclab.QRotation( theta4 );
      rot5 = qclab.QRotation( theta5 );
      rot6 = qclab.QRotation( theta6 );
      
      G5 = TFXY( [4,2], rot1, rot2, rot3, rot4, rot5, rot6 );
      test.verifyEqual( G5.matrix, G3.matrix, 'AbsTol', tol );
      
      % 6.
      rot = [rot1, rot2, rot3, rot4, rot5, rot6];
      G6 = TFXY( [4,2], rot );
      test.verifyTrue( G5 == G6 );
      
      % 7.
      TFXYMat = @f3c.qgates.RotationTFXYMatrix ;
      GM7 = TFXYMat( G6 );
      G7 = TFXY( GM7 );
      
      qubits = G7.qubits;
      test.verifyEqual( length(qubits), 2 );
      test.verifyEqual( qubits(1), int32(2) );
      test.verifyEqual( qubits(2), int32(4) );
      
      test.verifyEqual( G7.matrix, G6.matrix, 'AbsTol', tol );
      test.verifyEqual( G7.matrix, GM7.matrix, 'AbsTol', tol );
    end
    
    function test_RotationTFXY_update(test)
      tol = 10*eps;
      TFXY = @f3c.qgates.RotationTFXY ;
      
      G1 = TFXY() ;
      % update and angles check
      
      % 1.
      theta1 = pi/4;
      theta2 = -pi/4;
      theta3 = 5;
      theta4 = 0.1;
      theta5 = 5*pi/2;
      theta6 = pi/8;
      
      G1.update( theta1, theta2, theta3, theta4, theta5, theta6 );
      
      % left circuit angles
      lTheta = G1.thetasLeft ;
      modDif = mod(lTheta(1) - pi/4, 2*pi);
      test.verifyTrue( min(modDif, modDif - 2*pi) < tol  );
      modDif = mod(lTheta(2) + pi/4, 2*pi);
      test.verifyTrue( min(modDif, modDif - 2*pi) < tol  );
      % middle circuit angles
      mTheta = G1.thetasMiddle ;
      modDif = mod(mTheta(1) - 5, 2*pi);
      test.verifyTrue( min(modDif, modDif - 2*pi) < tol  );
      modDif = mod(mTheta(2) - 0.1, 2*pi);
      test.verifyTrue( min(modDif, modDif - 2*pi) < tol  );
      % right circuit angles
      rTheta = G1.thetasRight ;
      modDif = mod(rTheta(1) - 5*pi/2, 2*pi);
      test.verifyTrue( min(modDif, modDif - 2*pi) < tol  );
      modDif = mod(rTheta(2) - pi/8, 2*pi);
      test.verifyTrue( min(modDif, modDif - 2*pi) < tol  );
      
      % 2. 
      theta = [ theta1, theta2, theta3, theta4, theta5, theta6 ];
      G2 = TFXY();
      
      G2.update( theta );
      test.verifyTrue( G1 == G2 );
      
      % 3.
      rot1 = qclab.QRotation( theta1 );
      rot2 = qclab.QRotation( theta2 );
      rot3 = qclab.QRotation( theta3 );
      rot4 = qclab.QRotation( theta4 );
      rot5 = qclab.QRotation( theta5 );
      rot6 = qclab.QRotation( theta6 );
      
      G3 = TFXY() ;
      G3.update( rot1, rot2, rot3, rot4, rot5, rot6 );
      test.verifyEqual( G3.matrix, G2.matrix, 'AbsTol', tol );
      
      % 4.
      rot = [rot1, rot2, rot3, rot4, rot5, rot6];
      G4 = TFXY( );
      G4.update( rot )
      test.verifyTrue( G3 == G4 );
      
      % 5.
      TFXYMat = @f3c.qgates.RotationTFXYMatrix ;
      GM5 = TFXYMat( G4 );
      G5 = TFXY( );
      G5.update( GM5 );
      
      test.verifyEqual( G5.matrix, G4.matrix, 'AbsTol', tol );
      test.verifyEqual( G5.matrix, GM5.matrix, 'AbsTol', tol );
      
    end
    
    
    function test_RotationTFXY_fuse(test)
      
      tol = 10*eps;
      TFXY = @f3c.qgates.RotationTFXY ;
      
      G1 = TFXY( [0,1], pi/4, -pi/4, 5, 0.1, 5*pi/2, pi/8 ) ;
      G2 = TFXY( [0,1], 8*pi/4, pi/4, -7, pi, 9*pi/2, pi/8 ) ;
      G1c = copy(G1);
      G1.fuse( G2 );
      test.verifyEqual( G1.matrix, G2.matrix * G1c.matrix, 'AbsTol', tol );

    end
    
    function test_RotationTFXY_turnover( test )
      
      tol = 100*eps;
      TFXY = @f3c.qgates.RotationTFXY ;
      
      G1 = TFXY( [0,1], pi * randn(6,1) );
      G2 = TFXY( [1,2], pi * randn(6,1) );
      G3 = TFXY( [0,1], pi * randn(6,1) );
      [GA, GB, GC] = turnover( G1, G2, G3 );
      
      Cin = qclab.QCircuit( 3 );
      Cin.push_back( G1 );
      Cin.push_back( G2 );
      Cin.push_back( G3 );
      
      Cout = qclab.QCircuit( 3 );
      Cout.push_back( GA );
      Cout.push_back( GB );
      Cout.push_back( GC );
      
      test.verifyEqual( GA.qubit, int32(1) );
      test.verifyEqual( GB.qubit, int32(0) );
      test.verifyEqual( GC.qubit, int32(1) );
      
      test.verifyEqual(Cin.matrix, Cout.matrix, 'AbsTol', tol );
    end
    
    function test_RotationTFXY_draw_qasm( test )
      TFXY = @f3c.qgates.RotationTFXY ;
      G = TFXY( [0,1], pi * randn(6,1) );
      G.draw(1, 'N');
      G.draw(1, 'S');
      G.draw(1, 'L');
      
      G = TFXY( [4,5], 0.1, 0.2, 0.3, 0.4, 0.5, 0.6 );
      [T,out] = evalc('G.toQASM(1)'); % capture output to std::out in T
      test.verifyEqual( out, 0 );
      QASMstring = [sprintf('rz(%.15f) q[%d];\n', 0.1, 4), ...
                    sprintf('rz(%.15f) q[%d];\n', 0.2, 5 ), ...
                    sprintf('rx(%.15f) q[%d];\n', pi/2, 4 ), ...
                    sprintf('rx(%.15f) q[%d];\n', pi/2, 5 ), ...
                    sprintf('cx q[%d], q[%d];\n', 4, 5), ...
                    sprintf('rx(%.15f) q[%d];\n', 0.3, 4 ), ...
                    sprintf('rz(%.15f) q[%d];\n', 0.4, 5 ), ...
                    sprintf('cx q[%d], q[%d];\n', 4, 5), ...
                    sprintf('rx(%.15f) q[%d];\n', -pi/2, 4 ), ...
                    sprintf('rx(%.15f) q[%d];\n', -pi/2, 5 ), ...
                    sprintf('rz(%.15f) q[%d];\n', 0.5, 4), ...
                    sprintf('rz(%.15f) q[%d];\n', 0.6, 5 )];
      test.verifyEqual(T(1:length(QASMstring)), QASMstring);
    end
  end
end