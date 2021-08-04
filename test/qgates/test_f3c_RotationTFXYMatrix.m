classdef test_f3c_RotationTFXYMatrix < matlab.unittest.TestCase
  methods (Test)
    function test_RotationTFXYMatrix_basic_functionality(test)
      
      tol = 10*eps;
      TFXY = @f3c.qgates.RotationTFXY ;
      TFXYMat = @f3c.qgates.RotationTFXYMatrix ;
      
      % default constructor
      GM1 = TFXYMat();
      test.verifyEqual( GM1.nbQubits, int32(2) );    % nbQubits
      test.verifyFalse( GM1.fixed );                 % fixed
      test.verifyFalse( GM1.controlled );            % controlled
      
      % matrix
      test.verifyEqual( GM1.matrix, eye(4) );
      
      % qubits
      qubits = GM1.qubits;
      test.verifyEqual( length(qubits), 2 );
      test.verifyEqual( qubits(1), int32(0) );
      test.verifyEqual( qubits(2), int32(1) );
      qnew = [5, 3] ;
      GM1.setQubits( qnew );
      qubits = GM1.qubits;
      test.verifyEqual( qubits(1), int32(3) );
      test.verifyEqual( qubits(2), int32(5) );
      
      % fixed
      GM1.makeFixed();
      test.verifyTrue( GM1.fixed );
      GM1.makeVariable();
      test.verifyFalse( GM1.fixed );
      
      % only qubit constructor
      GM1 = TFXYMat( [5, 7] );
      test.verifyEqual( GM1.nbQubits, int32(2) );    % nbQubits
      qubits = GM1.qubits;
      test.verifyEqual( length(qubits), 2 );
      test.verifyEqual( qubits(1), int32(5) );
      test.verifyEqual( qubits(2), int32(7) );

      % qubit and vector constructor
      c1 = cos( pi/5) * exp(1i * 0.33);
      s1 = sin( pi/5) * exp(1i * 0.33);
      c2 = cos( 1.8994 ) * exp(-1i * pi/7);
      s2 = sin( 1.8994 ) * exp(-1i * pi/7);
      GM1 = TFXYMat( [1, 2], [c1, c2, s2, s1] );
      mat = GM1.matrix ;
      test.verifyEqual( mat(1,1), c1 );
      test.verifyEqual( mat(2,2), c2 );
      test.verifyEqual( mat(3,2), s2 );
      test.verifyEqual( mat(4,1), s1 );
      % cast from TFXYMat to TFXY gate
      G1 = TFXY( GM1 );
      test.verifyEqual( G1.matrix, GM1.matrix, 'AbsTol', tol );

      % cast from TFXY gate to TFXYMat
      G1 = TFXY([5, 3], -2*pi, -2*pi/2, -2*pi/3, 2*pi/5, 2*pi/4, 2*pi/3);
      GM1 = TFXYMat( G1 );
      
      % qubits
      qubits = GM1.qubits;
      test.verifyEqual( length(qubits), 2 );
      test.verifyEqual( qubits(1), int32(3) );
      test.verifyEqual( qubits(2), int32(5) );
      
      % matrix
      test.verifyEqual( G1.matrix, GM1.matrix, 'AbsTol', tol );
      
    end
    
    function test_RotationTFXYMatrix_fuse(test)
      
      tol = 10*eps;
      TFXY = @f3c.qgates.RotationTFXY ;
      TFXYMat = @f3c.qgates.RotationTFXYMatrix ;
      
      G1 = TFXY([0, 1], -2*pi, -2*pi/2, -2*pi/3, 2*pi/5, 2*pi/4, 2*pi/3);
      GM1 = TFXYMat( G1 );
      
      G2 = copy( G1 );
      GM2 = TFXYMat( G2 );
      
      G1c = copy(G1);
      
      % right fuse
      G1.fuse( G2, 'R' );
      test.verifyEqual( G1.matrix, G2.matrix * G1c.matrix, 'AbsTol', tol );
      
      GM1c = copy(GM1);
      GM1.fuse( GM2, 'R' );
      test.verifyEqual( GM1.matrix, GM2.matrix * GM1c.matrix, 'AbsTol', tol );
      
      test.verifyEqual( GM1.matrix, G1.matrix, 'AbsTol', tol );
      
      % left fuse
      G1 = copy(G1c);
      G1.fuse( G2, 'L' );
      test.verifyEqual( G1.matrix, G1c.matrix * G2.matrix, 'AbsTol', tol );
      
      GM1 = copy(GM1c);
      GM1.fuse( GM2, 'L' );
      test.verifyEqual( GM1.matrix, GM1c.matrix * GM2.matrix, 'AbsTol', tol );
      
      test.verifyEqual( GM1.matrix, G1.matrix, 'AbsTol', tol );
      
    end
    
    function test_RotationTFXYMatrix_turnover( test )
       tol = 10*eps;
       
       TFXY = @f3c.qgates.RotationTFXY ;
       TFXYMat = @f3c.qgates.RotationTFXYMatrix ;
       Circuit = @qclab.QCircuit ;
       
       % test Vee to Hat direction
       G1 = TFXY([0, 1], -2*pi, -2*pi/2, -2*pi/3, 2*pi/5, 2*pi/4, 2*pi/3);
       GM1 = TFXYMat( G1 );
       
       G2 = TFXY([1, 2], 0.1, 0.2, -0.23, 0.45, 1.2, -3.11) ;
       GM2 = TFXYMat( G2 );
       
       G3 = TFXY([0, 1], -0.45, -2.43, 1.11, 2.76, -2.94, -1.17) ;
       GM3 = TFXYMat( G3 );
       
       % check qubits
       qubits1 = GM1.qubits;
       test.verifyEqual( length(qubits1), 2 );
       test.verifyEqual( qubits1(1), int32(0) );
       test.verifyEqual( qubits1(2), int32(1) );
       
       qubits2 = GM2.qubits;
       test.verifyEqual( length(qubits2), 2 );
       test.verifyEqual( qubits2(1), int32(1) );
       test.verifyEqual( qubits2(2), int32(2) );
       
       qubits3 = GM3.qubits;
       test.verifyEqual( length(qubits3), 2 );
       test.verifyEqual( qubits3(1), int32(0) );
       test.verifyEqual( qubits3(2), int32(1) );
       
       % check circuit matrices before turnover
       C_Theta = Circuit( 3 ); 
       C_Theta.push_back( G1 );
       C_Theta.push_back( G2 );
       C_Theta.push_back( G3 );
       
       C_Mat = Circuit( 3 );
       C_Mat.push_back( GM1 );
       C_Mat.push_back( GM2 );
       C_Mat.push_back( GM3 );
       
       Mat = C_Mat.matrix ;
       
       test.verifyEqual( C_Theta.matrix, Mat, 'AbsTol', tol ) ;
       
       % evaluate turnover operation on pattern
       [GMA, GMB, GMC] = turnover( GM1, GM2, GM3 );

       C_Mat.replace( 1, GMA );
       C_Mat.replace( 2, GMB );
       C_Mat.replace( 3, GMC );

       % check qubits after turnover
       qubitsA = GMA.qubits;
       test.verifyEqual( qubitsA, qubits2 );
       qubitsB = GMB.qubits;
       test.verifyEqual( qubitsB, qubits1 );  
       qubitsC = GMC.qubits;
       test.verifyEqual( qubitsC, qubits2 );
       
       % check circuit after turnover and compare with before
       test.verifyEqual( C_Mat.matrix, Mat, 'AbsTol', tol ) ;
       
       % test Hat to Vee direction
       G1 = TFXY([1, 2], -2*pi, -2*pi/2, -2*pi/3, 2*pi/5, 2*pi/4, 2*pi/3);
       GM1 = TFXYMat( G1 );
       
       G2 = TFXY([0, 1], 0.1, 0.2, -0.23, 0.45, 1.2, -3.11) ;
       GM2 = TFXYMat( G2 );
       
       G3 = TFXY([1, 2], -0.45, -2.43, 1.11, 2.76, -2.94, -1.17) ;
       GM3 = TFXYMat( G3 );
       
       % check qubits
       qubits1 = GM1.qubits;
       test.verifyEqual( length(qubits1), 2 );
       test.verifyEqual( qubits1(1), int32(1) );
       test.verifyEqual( qubits1(2), int32(2) );
       
       qubits2 = GM2.qubits;
       test.verifyEqual( length(qubits2), 2 );
       test.verifyEqual( qubits2(1), int32(0) );
       test.verifyEqual( qubits2(2), int32(1) );
       
       qubits3 = GM3.qubits;
       test.verifyEqual( length(qubits3), 2 );
       test.verifyEqual( qubits3(1), int32(1) );
       test.verifyEqual( qubits3(2), int32(2) );
       
       % check circuit matrices before turnover
       C_Theta = Circuit( 3 ); 
       C_Theta.push_back( G1 );
       C_Theta.push_back( G2 );
       C_Theta.push_back( G3 );
       
       C_Mat = Circuit( 3 );
       C_Mat.push_back( GM1 );
       C_Mat.push_back( GM2 );
       C_Mat.push_back( GM3 );
       
       Mat = C_Mat.matrix ;
       
       test.verifyEqual( C_Theta.matrix, Mat, 'AbsTol', tol ) ;
       
       % evaluate turnover operation on pattern
       C_Mat.draw( 1, 'S' );
       [GMA, GMB, GMC] = turnover( GM1, GM2, GM3 );
       C_Mat.replace( 1, GMA );
       C_Mat.replace( 2, GMB );
       C_Mat.replace( 3, GMC );
       C_Mat.draw( 1, 'S' );
       
       % check qubits after turnover
       qubitsA = GMA.qubits;
       test.verifyEqual( qubitsA, qubits2 );
       qubitsB = GMB.qubits;
       test.verifyEqual( qubitsB, qubits1 );  
       qubitsC = GMC.qubits;
       test.verifyEqual( qubitsC, qubits2 );
       
       % check circuit after turnover and compare with before
       test.verifyEqual( C_Mat.matrix, Mat, 'AbsTol', tol ) ;
    end
    
    
    function test_RotationTFXYMatrix_turnover_edge_cases( test )
      tol = 100 * eps ;
      test.runEdgeTurnoverTest( [1e-3, 1e3], true, tol, test );
      test.runEdgeTurnoverTest( [1e-3, 1e3], false, tol, test );
      
      test.runEdgeTurnoverTest( [1e-7, 1e7], true, tol, test );
      test.runEdgeTurnoverTest( [1e-7, 1e7], false, tol, test );
      
      test.runEdgeTurnoverTest( [1e-12, 1e12], true, tol, test );
      test.runEdgeTurnoverTest( [1e-12, 1e12], false, tol, test );
      
      test.runEdgeTurnoverTest( [1e-1, 1e3], true, tol, test );
      test.runEdgeTurnoverTest( [1e-1, 1e3], false, tol, test );
      
      test.runEdgeTurnoverTest( [1e-3, 1e1], true, tol, test );
      test.runEdgeTurnoverTest( [1e-3, 1e1], false, tol, test );
    end
    
    function test_RotationTFXYMatrix_turnover_difficult_cases( test )
      tol = 100 * eps ;
      
      % ------------------------------------------------------------------------
      abcd1 = [9.977474154852748e-01 + 6.708274660786033e-02i, ...     
         0.000000000000000e+00 + 0.000000000000000e+00i, ...     
         5.747921180094008e-01 - 8.182994690663481e-01i, ... 
         6.028190300314134e-10 - 7.978779461997844e-10i];
       
       abcd2 = [0.000000000000000e+00 + 0.000000000000000e+00i, ...
         7.188471818119035e-01 - 6.951681301678638e-01i, ...
         2.478843595396635e-10 + 9.687896285033252e-10i, ...
         9.577767742158577e-01 - 2.875128706208922e-01i];
       
       abcd3 = [0.000000000000000e+00 + 0.000000000000000e+00i, ...
         0.000000000000000e+00 + 0.000000000000000e+00i, ...
         9.581348980296205e-01 + 2.863171618638478e-01i, ...
         9.999599996633471e-01 + 8.944220104563980e-03i];
       
      gate1 = f3c.qgates.RotationTFXYMatrix( [0,1], abcd1 );
      gate2 = f3c.qgates.RotationTFXYMatrix( [1,2], abcd2 );
      gate3 = f3c.qgates.RotationTFXYMatrix( [0,1], abcd3 );
      test.runTurnoverTest( gate1, gate2, gate3, tol, test ) ;
      % ------------------------------------------------------------------------
      
      abcd1 = [0.000000000000000e+00 + 0.000000000000000e+00i, ...
         0.000000000000000e+00 + 0.000000000000000e+00i, ...
         6.647942024722138e-01 - 7.470265513081401e-01i, ...
         8.461608114044463e-01 + 5.329276510403349e-01i];
       
      abcd2 = [0.000000000000000e+00 + 0.000000000000000e+00i, ...
         7.506277343123036e-01 - 6.607253623716725e-01i, ...
         9.736658994431211e-10 + 2.279796400155463e-10i, ...
         8.582907146965845e-01 + 5.131637643731544e-01i ];
       
      abcd3 = [9.637370846311818e-01 + 2.668535772789833e-01i, ...
         0.000000000000000e+00 + 0.000000000000000e+00i, ...
         -2.882863951963193e-01 - 9.575442310121823e-01i, ...
         2.417758192331365e-10 + 9.703321355258445e-10i ];
       
      gate1 = f3c.qgates.RotationTFXYMatrix( [0,1], abcd1 );
      gate2 = f3c.qgates.RotationTFXYMatrix( [1,2], abcd2 );
      gate3 = f3c.qgates.RotationTFXYMatrix( [0,1], abcd3 );
      test.runTurnoverTest( gate1, gate2, gate3, tol, test ) ;
      % ------------------------------------------------------------------------
      
      abcd1 = [0.000000000000000e+00 + 0.000000000000000e+00i, ...
         0.000000000000000e+00 + 0.000000000000000e+00i, ...
         8.722079837813873e-01 - 4.891351889079412e-01i, ...
         9.114066170405436e-01 + 4.115069603478316e-01i];

  
       abcd2 = [0.000000000000000e+00 + 0.000000000000000e+00i, ...
         0.000000000000000e+00 + 0.000000000000000e+00i, ...
         9.366150910289665e-01 - 3.503600594485634e-01i, ...
         1.510323145929456e-01 - 9.885288260585513e-01i];

       abcd3 = [ 0.000000000000000e+00 + 0.000000000000000e+00i, ...
          5.228005520286667e-01 - 8.524550327134689e-01i, ...
          4.230742476926466e-10 + 9.060950176109023e-10i, ...
          7.608367182507818e-01 - 6.489433628302090e-01i ];
        
      gate1 = f3c.qgates.RotationTFXYMatrix( [0,1], abcd1 );
      gate2 = f3c.qgates.RotationTFXYMatrix( [1,2], abcd2 );
      gate3 = f3c.qgates.RotationTFXYMatrix( [0,1], abcd3 );
      test.runTurnoverTest( gate1, gate2, gate3, tol, test ) ;
      % ------------------------------------------------------------------------
      
      abcd1 = [
      0.000000000000000e+00 + 0.000000000000000e+00i
      7.565770258663967e-01 - 6.539045832009114e-01i
      8.599866250772952e-13 + 5.103165730095034e-13i
      9.959933714808024e-01 + 8.942708743051155e-02i
      ];
    
      abcd2 = [
      6.078637054323632e-01 + 7.940413815526475e-01i
      9.325857931676858e-01 - 3.609483874209697e-01i
      3.122332486560884e-13 - 9.500054728440595e-13i
      6.663868822784031e-13 - 7.456061447756916e-13i
      ];
    
      abcd3 = [
      0.000000000000000e+00 + 0.000000000000000e+00i
      8.111316659087202e-01 + 5.848635914126852e-01i
      9.363477918049556e-13 - 3.510737996233608e-13i
     -1.166118067038798e-01 + 9.931775704964632e-01i
      ];
    
      gate1 = f3c.qgates.RotationTFXYMatrix( [0,1], abcd1 );
      gate2 = f3c.qgates.RotationTFXYMatrix( [1,2], abcd2 );
      gate3 = f3c.qgates.RotationTFXYMatrix( [0,1], abcd3 );
      test.runTurnoverTest( gate1, gate2, gate3, tol, test ) ;
      % ------------------------------------------------------------------------
      
      abcd1 = [
      -7.414033686368241e-03 - 9.999725156745546e-01i
      0.000000000000000e+00 + 0.000000000000000e+00i
      9.779047952095287e-01 - 2.090507390711878e-01i
      6.339403396197329e-13 - 7.733819533728580e-13i
      ];
    
      abcd2 = [
      7.383926033647923e-01 - 6.743710872332567e-01i
      1.303603212772043e-01 + 9.914666845822425e-01i
      9.993576134117824e-13 - 3.583797589020022e-14i
      6.898383132047741e-13 - 7.239634670581051e-13i
      ];

      abcd3 = [
      -3.139884909362514e-01 - 9.494267889413989e-01i
      0.000000000000000e+00 + 0.000000000000000e+00i
      5.918465684454635e-01 + 8.060506432100463e-01i
      4.140852982486016e-13 - 9.102380819183333e-13i
      ];
    
      gate1 = f3c.qgates.RotationTFXYMatrix( [0,1], abcd1 );
      gate2 = f3c.qgates.RotationTFXYMatrix( [1,2], abcd2 );
      gate3 = f3c.qgates.RotationTFXYMatrix( [0,1], abcd3 );
      test.runTurnoverTest( gate1, gate2, gate3, tol, test ) ;
      % ------------------------------------------------------------------------
      
      abcd1 = [
      9.769835341012577e-01 - 2.133147301407398e-01i
      3.167717394825452e-01 + 9.485018002435222e-01i
      6.810075261366204e-13 - 7.322764159422863e-13i
      4.404757727317805e-13 + 8.977644978703161e-13i
      ];

      abcd2 = [
      0.000000000000000e+00 + 0.000000000000000e+00i
      8.070981373579820e-01 + 5.904173072270799e-01i
     -5.825078059787758e-13 - 8.128251078637967e-13i
      9.943275242318647e-01 + 1.063615276071688e-01i
      ];
   
      abcd3 = [
      8.516832738839737e-01 - 5.240568680842531e-01i
      7.013626296997841e-01 + 7.128046448085221e-01i
      8.968925144185521e-13 + 4.422485925132689e-13i
     -3.386881907940717e-13 + 9.408986711738084e-13i
      ];
    
      gate1 = f3c.qgates.RotationTFXYMatrix( [0,1], abcd1 );
      gate2 = f3c.qgates.RotationTFXYMatrix( [1,2], abcd2 );
      gate3 = f3c.qgates.RotationTFXYMatrix( [0,1], abcd3 );
      test.runTurnoverTest( gate1, gate2, gate3, tol, test ) ;
      % ------------------------------------------------------------------------
      
      abcd1 = [
      9.625118278745912e-01 - 2.712397116970768e-01i
      5.600756247208285e-01 - 8.284414853165996e-01i
      8.327329076804433e-13 + 5.536749086477319e-13i
     -5.867691529444130e-14 + 9.982770254851751e-13i
      ];

      abcd2 = [
      7.647726157454231e-01 + 6.443002764285480e-01i
      0.000000000000000e+00 + 0.000000000000000e+00i
      5.089390435208107e-01 + 8.608025615552747e-01i
      9.699842591682075e-13 - 2.431677136585032e-13i
      ];

      abcd3 = [
      4.656330623995300e-01 + 8.849778817577507e-01i
      9.099900733823901e-01 - 4.146300355081771e-01i
     -1.952845598202712e-13 + 9.807466241062484e-13i
      6.404274666559924e-13 + 7.680186585966435e-13i
      ];
    
      gate1 = f3c.qgates.RotationTFXYMatrix( [0,1], abcd1 );
      gate2 = f3c.qgates.RotationTFXYMatrix( [1,2], abcd2 );
      gate3 = f3c.qgates.RotationTFXYMatrix( [0,1], abcd3 );
      test.runTurnoverTest( gate1, gate2, gate3, tol, test ) ;
      % ------------------------------------------------------------------------
      
      hat = false;
      a1 =  7.453921811940983e-08 + 1i * -6.491430950759949e-08;
      b1 =  9.556230583231934e-08 + 1i *  2.525486726774946e-08;
      c1 = -4.834527059706285e-02 + 1i * -9.988306837552033e-01;
      d1 =  5.697046993173709e-01 + 1i *  8.218494725773657e-01;
      a2 =  7.830553348877711e-08 + 1i * -6.031754339940732e-08;
      b2 =  9.550976184196558e-01 + 1i *  2.962913081565026e-01;
      c2 =  8.812551919021028e-08 + 1i * -4.726407586588123e-08;
      d2 =  3.078279795115421e-01 + 1i * -9.514420292534022e-01;
      a3 =  8.279991835076616e-08 + 1i * -5.398273925809837e-08;
      b3 =  9.804799194517483e-08 + 1i *  1.251214978412220e-08;
      c3 =  9.872727849708168e-01 + 1i * -1.590359961013784e-01;
      d3 =  9.203567108822643e-01 + 1i *  3.910799467295281e-01;
      
      if hat 
        q13 = [1,2];
        q2 = [0,1];
      else
        q13 = [0,1];
        q2 = [1,2];
      end
      
      gate1 = f3c.qgates.RotationTFXYMatrix( q13, [a1, b1, c1, d1] );
      gate2 = f3c.qgates.RotationTFXYMatrix( q2, [a2, b2, c2, d2] );
      gate3 = f3c.qgates.RotationTFXYMatrix( q13, [a3, b3, c3, d3] );
      test.runTurnoverTest( gate1, gate2, gate3, tol, test ) ;
      % ------------------------------------------------------------------------
      
      hat = true;
      a1 =  6.308193058238803e-01 + 1i *  7.759297670536088e-01;
      b1 = -1.226359552239620e-01 + 1i *  9.924517230003161e-01;
      c1 =  8.189600410805670e-08 + 1i *  5.738505477154360e-08;
      d1 =  2.260477550526254e-08 + 1i *  9.741162212157532e-08;
      a2 =  9.938184330411782e-01 + 1i * -1.110176659255953e-01;
      b2 =  9.020255187189186e-08 + 1i *  4.041611377285099e-08;
      c2 =  8.814403320081540e-01 + 1i * -4.722953960281059e-01;
      d2 =  4.023943241383570e-08 + 1i * -9.154664428045545e-08;
      a3 =  8.995987376780848e-01 + 1i *  4.367174271402349e-01;
      b3 =  8.208849271282048e-01 + 1i * -5.710936319148655e-01;
      c3 =  9.736967705038767e-08 + 1i *  2.278477542358291e-08;
      d3 =  7.468406901538925e-08 + 1i *  6.650029951289297e-08;
      
      if hat 
        q13 = [1,2];
        q2 = [0,1];
      else
        q13 = [0,1];
        q2 = [1,2];
      end
      
      gate1 = f3c.qgates.RotationTFXYMatrix( q13, [a1, b1, c1, d1] );
      gate2 = f3c.qgates.RotationTFXYMatrix( q2, [a2, b2, c2, d2] );
      gate3 = f3c.qgates.RotationTFXYMatrix( q13, [a3, b3, c3, d3] );
      test.runTurnoverTest( gate1, gate2, gate3, tol, test ) ;
      % ------------------------------------------------------------------------
      
      hat = false;
      a1 =  8.040304886458259e-08 + 1i * -5.749184594341116e-08;
      b1 =  7.373138834590875e-08 + 1i *  6.583042601477045e-08;
      c1 =  9.983858291235649e-01 + 1i * -5.679556501384782e-02;
      d1 =  9.985685513389021e-01 + 1i * -5.348689817999129e-02;
      a2 = -5.169668546804159e-08 + 1i *  8.424615913096427e-08;
      b2 =  9.595972106656034e-01 + 1i *  2.813773147977348e-01;
      c2 =  4.051146504038741e-08 + 1i * -9.142658913183499e-08;
      d2 =  2.800326978703757e-01 + 1i *  9.599904625169091e-01;
      a3 =  7.450242420396489e-01 + 1i * -6.670373893367871e-01;
      b3 =  1.602724700240165e-01 + 1i *  9.870728115759195e-01;
      c3 =  5.914974679167009e-08 + 1i * -8.063068556375551e-08;
      d3 =  9.992877767346097e-08 + 1i *  3.773511983273288e-09;
      
      if hat
        q13 = [1,2];
        q2 = [0,1];
      else
        q13 = [0,1];
        q2 = [1,2];
      end
      
      gate1 = f3c.qgates.RotationTFXYMatrix( q13, [a1, b1, c1, d1] );
      gate2 = f3c.qgates.RotationTFXYMatrix( q2, [a2, b2, c2, d2] );
      gate3 = f3c.qgates.RotationTFXYMatrix( q13, [a3, b3, c3, d3] );
      test.runTurnoverTest( gate1, gate2, gate3, tol, test ) ;
      % ------------------------------------------------------------------------
      
      hat = true;
      a1 =  8.508996043460090e-01 + 1i * -5.253283385881590e-01;
      b1 =  9.524673879925145e-01 + 1i * -3.046405665874241e-01;
      c1 =  9.822975476171381e-08 + 1i * -1.873273283462558e-08;
      d1 =  9.995811443211189e-08 + 1i * -2.894021281318401e-09;
      a2 =  9.750362944107126e-01 + 1i * -2.220455461877498e-01;
      b2 =  3.462360770782283e-01 + 1i * -9.381474185486358e-01;
      c2 =  1.887027753396164e-08 + 1i * -9.820342471518578e-08;
      d2 =  5.175058712180079e-08 + 1i *  8.556796557444150e-08;
      a3 =  7.057856552777438e-08 + 1i *  6.919991838682351e-08;
      b3 = -9.433603801253276e-09 + 1i *  9.839191905854035e-08;
      c3 =  9.994792221685084e-01 + 1i * -3.226893945303590e-02;
      d3 =  3.044015169411465e-01 + 1i *  9.525438134195818e-01;

      if hat
        q13 = [1,2];
        q2 = [0,1];
      else
        q13 = [0,1];
        q2 = [1,2];
      end
      
      gate1 = f3c.qgates.RotationTFXYMatrix( q13, [a1, b1, c1, d1] );
      gate2 = f3c.qgates.RotationTFXYMatrix( q2, [a2, b2, c2, d2] );
      gate3 = f3c.qgates.RotationTFXYMatrix( q13, [a3, b3, c3, d3] );
      test.runTurnoverTest( gate1, gate2, gate3, tol, test ) ;
      % ------------------------------------------------------------------------ 
      
      hat = true ;
      a1 =  6.564708629609297e-08 + 1i * -7.389467286303323e-08;
      b1 =  3.243734243525461e-08 + 1i *  9.336906035962560e-08;
      c1 =  9.361911630527955e-01 + 1i *  3.514912605198661e-01;
      d1 =  9.885820908071820e-01 + 1i *  1.506832762296171e-01;
      a2 =  8.705659527048012e-08 + 1i *  4.680931356699430e-08;
      b2 =  3.396895062872275e-01 + 1i * -9.405376331217852e-01;
      c2 =  9.523145319934191e-08 + 1i * -3.051180626481246e-08;
      d2 =  9.602999168943377e-01 + 1i *  2.789696571541758e-01;
      a3 =  3.468261723642508e-08 + 1i *  9.255851488832893e-08;
      b3 =  7.126729740220410e-01 + 1i *  7.014964234396133e-01;
      c3 =  9.956483366461913e-08 + 1i * -9.319007314978578e-09;
      d3 =  1.820095624442049e-01 + 1i * -9.832967604842697e-01;
      
      if hat
        q13 = [1,2];
        q2 = [0,1];
      else
        q13 = [0,1];
        q2 = [1,2];
      end
      
      gate1 = f3c.qgates.RotationTFXYMatrix( q13, [a1, b1, c1, d1] );
      gate2 = f3c.qgates.RotationTFXYMatrix( q2, [a2, b2, c2, d2] );
      gate3 = f3c.qgates.RotationTFXYMatrix( q13, [a3, b3, c3, d3] );
      test.runTurnoverTest( gate1, gate2, gate3, tol, test ) ;
      % ------------------------------------------------------------------------ 
      
      hat = true;
      a1 = 6.340924154984683e-01 + 1i * -7.732572719395001e-01;
      b1 = 0.000000000000000e+00 + 1i * 0.000000000000000e+00;
      c1 = 9.713077446841060e-01 + 1i * -2.378261237052726e-01;
      d1 = 9.999532192647693e-13 + 1i * 9.672604717653607e-15;
      a2 = 0.000000000000000e+00 + 1i * 0.000000000000000e+00;
      b2 = 0.000000000000000e+00 + 1i * 0.000000000000000e+00;
      c2 = 9.471195219713482e-01 + 1i * 3.208809921150908e-01;
      d2 = 8.690578371831881e-01 + 1i * 4.947104967862310e-01;
      a3 = 0.000000000000000e+00 + 1i * -0.000000000000000e+00;
      b3 = 8.761394757430743e-01 + 1i * -4.820576926516690e-01;
      c3 = 9.788529481287986e-13 + 1i * -2.045651630643877e-13;
      d3 = 9.994141846336670e-01 + 1i * -3.422407855622338e-02;
      
      if hat
        q13 = [1,2];
        q2 = [0,1];
      else
        q13 = [0,1];
        q2 = [1,2];
      end
      
      gate1 = f3c.qgates.RotationTFXYMatrix( q13, [a1, b1, c1, d1] );
      gate2 = f3c.qgates.RotationTFXYMatrix( q2, [a2, b2, c2, d2] );
      gate3 = f3c.qgates.RotationTFXYMatrix( q13, [a3, b3, c3, d3] );
      test.runTurnoverTest( gate1, gate2, gate3, tol, test ) ;
      % ------------------------------------------------------------------------ 
      
      hat = true;
      a1 =  5.597260981924069e-01 + 1i * -8.286776785954199e-01;
      b1 =  9.677577706900544e-01 + 1i *  2.518827053789244e-01;
      c1 =  9.478166376834505e-08 + 1i *  3.188159678096891e-08;
      d1 =  7.941370552610752e-08 + 1i * -6.077387082137010e-08;
      a2 =  7.278876416526158e-01 + 1i *  6.856964205312605e-01;
      b2 =  6.746369175766853e-02 + 1i *  9.977217298898607e-01;
      c2 =  2.761296116729983e-08 + 1i *  9.611204074190232e-08;
      d2 =  9.661397083863748e-08 + 1i *  2.580195029044911e-08;
      a3 = -5.034454048154690e-09 + 1i * -9.871482608574254e-08;
      b3 =  9.892381829124894e-01 + 1i *  1.463141054990464e-01;
      c3 =  4.327581657942693e-08 + 1i *  9.015100498266059e-08;
      d3 =  9.936461727682634e-01 + 1i * -1.125490263972667e-01;
      
      if hat
        q13 = [1,2];
        q2 = [0,1];
      else
        q13 = [0,1];
        q2 = [1,2];
      end
      
      gate1 = f3c.qgates.RotationTFXYMatrix( q13, [a1, b1, c1, d1] );
      gate2 = f3c.qgates.RotationTFXYMatrix( q2, [a2, b2, c2, d2] );
      gate3 = f3c.qgates.RotationTFXYMatrix( q13, [a3, b3, c3, d3] );
      test.runTurnoverTest( gate1, gate2, gate3, tol, test ) ;
      % ------------------------------------------------------------------------
      
      hat = true;
      a1 =  8.771666498303764e-01 + 1i * -4.801860768757814e-01;
      b1 =  0.000000000000000e+00 + 1i *  0.000000000000000e+00;
      c1 =  1.120006420710291e-01 + 1i *  9.937081343008506e-01;
      d1 =  5.033522803015200e-13 + 1i * -8.640812935802164e-13;
      a2 =  0.000000000000000e+00 + 1i *  0.000000000000000e+00;
      b2 =  0.000000000000000e+00 + 1i * -0.000000000000000e+00;
      c2 =  7.827050908579476e-01 + 1i * -6.223927544123983e-01;
      d2 =  8.697745471024860e-01 + 1i *  4.934493258812553e-01;
      a3 =  0.000000000000000e+00 + 1i *  0.000000000000000e+00;
      b3 =  9.982129849640763e-01 + 1i *  5.975647788406616e-02;
      c3 = -5.855328274366165e-14 + 1i *  9.982842847004760e-13;
      d3 = -2.079633343912144e-01 + 1i * -9.781366221284673e-01;
      
      if hat
        q13 = [1,2];
        q2 = [0,1];
      else
        q13 = [0,1];
        q2 = [1,2];
      end
      
      gate1 = f3c.qgates.RotationTFXYMatrix( q13, [a1, b1, c1, d1] );
      gate2 = f3c.qgates.RotationTFXYMatrix( q2, [a2, b2, c2, d2] );
      gate3 = f3c.qgates.RotationTFXYMatrix( q13, [a3, b3, c3, d3] );
      test.runTurnoverTest( gate1, gate2, gate3, tol, test ) ;
      % ------------------------------------------------------------------------
      hat = true;
      abcd1 = [
        0.000000000000000e+00 + 0.000000000000000e+00i
        0.000000000000000e+00 + 0.000000000000000e+00i
        -9.999875000260416e-01 - 4.999979166692711e-03i
        -2.499739591471234e-02 - 9.996875162757026e-01i
      ];
    
      abcd2 = [
        9.999500004166653e-01 + 3.469446951953614e-18i
        9.997500170828263e-01 - 1.999766676833116e-02i
        9.999833334166664e-03 + 0.000000000000000e+00i
        -9.997833434164497e-03 + 1.999833338388808e-04i

      ];
  
      abcd3 = [
        9.998708374009809e-03 - 1.499918751695294e-04i
        -4.999895834123260e-05 - 9.999708336510399e-03i
        9.999375010677004e-01 - 4.999729169817691e-03i
        1.499868754070252e-02 + 9.998375081508720e-01i
      ];
    
      if hat
        q13 = [1,2];
        q2 = [0,1];
      else
        q13 = [0,1];
        q2 = [1,2];
      end
      
      gate1 = f3c.qgates.RotationTFXYMatrix( q13, abcd1 );
      gate2 = f3c.qgates.RotationTFXYMatrix( q2, abcd2 );
      gate3 = f3c.qgates.RotationTFXYMatrix( q13, abcd3 );
      test.runTurnoverTest( gate1, gate2, gate3, tol, test ) ;
      % ------------------------------------------------------------------------
      hat = true;
    
      abcd1 = [
        -1.999867068204558e-02 + 9.998000065867929e-01i
        9.998000065867929e-01 - 1.999867068204558e-02i
        0.000000000000000e+00 + 0.000000000000000e+00i
        0.000000000000000e+00 + 0.000000000000000e+00i
      ];
    
      abcd2 = [
       9.997998681901258e-01 - 1.999866392343471e-02i
       9.999998614958482e-01 + 0.000000000000000e+00i
       5.263157651747099e-04 + 0.000000000000000e+00i
       -5.262105055303999e-04 + 1.052561356317559e-05i
      ];
  
      abcd3 = [
        -2.999550003570902e-02 - 9.995498951874905e-01i
        9.995498951874904e-01 + 2.999550003570902e-02i
        1.578710673967206e-05 - 5.260789407800226e-04i
        -5.260789407800226e-04 + 1.578710673967206e-05i
      ];
    
      if hat
        q13 = [1,2];
        q2 = [0,1];
      else
        q13 = [0,1];
        q2 = [1,2];
      end
      
      gate1 = f3c.qgates.RotationTFXYMatrix( q13, abcd1 );
      gate2 = f3c.qgates.RotationTFXYMatrix( q2, abcd2 );
      gate3 = f3c.qgates.RotationTFXYMatrix( q13, abcd3 );
      test.runTurnoverTest( gate1, gate2, gate3, tol, test ) ;
      % ------------------------------------------------------------------------
      hat = false;
    
      abcd1 = [
        9.997998681901258e-01 - 1.999866392343471e-02i
        9.999998614958482e-01 + 0.000000000000000e+00i
        5.263157651747099e-04 + 0.000000000000000e+00i
        -5.262105055303999e-04 + 1.052561356317559e-05i
      ];
    
      abcd2 = [
        -2.999550003570902e-02 - 9.995498951874905e-01i
        9.995498951874904e-01 + 2.999550003570902e-02i
        1.578710673967206e-05 - 5.260789407800226e-04i
        -5.260789407800226e-04 + 1.578710673967206e-05i
      ];
  
      abcd3 = [
        9.999494464277966e-01 - 9.999827794093308e-03i
        9.999994459834307e-01 + 0.000000000000000e+00i
        0.000000000000000e+00 - 1.052631384555583e-03i
        -1.052613840786892e-05 - 1.052578753424950e-03i
      ];
    
      if hat
        q13 = [1,2];
        q2 = [0,1];
      else
        q13 = [0,1];
        q2 = [1,2];
      end
      
      gate1 = f3c.qgates.RotationTFXYMatrix( q13, abcd1 );
      gate2 = f3c.qgates.RotationTFXYMatrix( q2, abcd2 );
      gate3 = f3c.qgates.RotationTFXYMatrix( q13, abcd3 );
      test.runTurnoverTest( gate1, gate2, gate3, tol, test ) ;
      % ------------------------------------------------------------------------
      hat = false;
      abcd1 = [
        9.999919996819380e-01 - 3.999989332026944e-03i
        9.999999996712687e-01 - 4.336808689942018e-19i
        2.564102563821597e-05 + 0.000000000000000e+00i
       -2.564082051028437e-05 + 1.025638290488092e-07i
      ];
    
      abcd2 = [
        9.999879641483611e-01 - 4.906210447323300e-03i
        4.906210447323300e-03 - 9.999879641483612e-01i
       -2.564071703506443e-05 - 1.258002679076543e-07i
       -1.258002679076543e-07 - 2.564071703506443e-05i
      ];
  
      abcd3 = [
        9.999979986857449e-01 - 1.999998664037086e-03i
        9.999999986850756e-01 + 0.000000000000000e+00i
        0.000000000000000e+00 - 5.128205125957394e-05i
        -1.025640341430932e-07 - 5.128194869550560e-05i
      ];
    
      if hat
        q13 = [1,2];
        q2 = [0,1];
      else
        q13 = [0,1];
        q2 = [1,2];
      end
      
      gate1 = f3c.qgates.RotationTFXYMatrix( q13, abcd1 );
      gate2 = f3c.qgates.RotationTFXYMatrix( q2, abcd2 );
      gate3 = f3c.qgates.RotationTFXYMatrix( q13, abcd3 );
      test.runTurnoverTest( gate1, gate2, gate3, tol, test ) ;
      % ------------------------------------------------------------------------
      hat = false;
      abcd1 = [
        9.999980000006666e-01 - 1.999998666666933e-03i
        1.000000000000000e+00 + 0.000000000000000e+00i
        0.000000000000000e+00 + 0.000000000000000e+00i
        0.000000000000000e+00 + 0.000000000000000e+00i
      ];
    
      abcd2 = [
        9.999995000000417e-01 - 9.999998333333417e-04i
        9.999995000000417e-01 + 9.999998333333417e-04i
        0.000000000000000e+00 + 0.000000000000000e+00i
        0.000000000000000e+00 + 0.000000000000000e+00i
      ];
  
      abcd3 = [
        9.999980000006666e-01 - 1.999998666666933e-03i
        9.999999997959392e-01 + 0.000000000000000e+00i
        0.000000000000000e+00 - 2.020202020064605e-05i
        0.000000000000000e+00 + 0.000000000000000e+00i
      ];
    
      if hat
        q13 = [1,2];
        q2 = [0,1];
      else
        q13 = [0,1];
        q2 = [1,2];
      end
      
      gate1 = f3c.qgates.RotationTFXYMatrix( q13, abcd1 );
      gate2 = f3c.qgates.RotationTFXYMatrix( q2, abcd2 );
      gate3 = f3c.qgates.RotationTFXYMatrix( q13, abcd3 );
      test.runTurnoverTest( gate1, gate2, gate3, tol, test ) ;
      % ------------------------------------------------------------------------
      hat = false;
      abcd1 = [
        9.950041652780257e-01 - 9.983341664682815e-02i
        1.000000000000000e+00 + 0.000000000000000e+00i
        0.000000000000000e+00 + 0.000000000000000e+00i
        0.000000000000000e+00 + 0.000000000000000e+00i
      ];
    
      abcd2 = [
        9.987502603949663e-01 - 4.997916927067833e-02i
        9.987502603949663e-01 + 4.997916927067833e-02i
        0.000000000000000e+00 + 0.000000000000000e+00i
        0.000000000000000e+00 + 0.000000000000000e+00i
      ];
  
      abcd3 = [
        9.950029997433866e-01 - 9.983329970329234e-02i
        9.999998698459004e-01 + 0.000000000000000e+00i
        0.000000000000000e+00 - 5.102040594976018e-04i
        -1.528061903045610e-04 - 1.522964964438442e-03i
      ];
    
      if hat
        q13 = [1,2];
        q2 = [0,1];
      else
        q13 = [0,1];
        q2 = [1,2];
      end
      
      gate1 = f3c.qgates.RotationTFXYMatrix( q13, abcd1 );
      gate2 = f3c.qgates.RotationTFXYMatrix( q2, abcd2 );
      gate3 = f3c.qgates.RotationTFXYMatrix( q13, abcd3 );
      test.runTurnoverTest( gate1, gate2, gate3, tol, test ) ;
      % ------------------------------------------------------------------------
      hat = false;
      
      abcd1 = [
        1.000000000000000e+00 + 0.000000000000000e+00i
        9.950041652780259e-01 - 9.983341664682817e-02i
        0.000000000000000e+00 + 0.000000000000000e+00i
        0.000000000000000e+00 + 0.000000000000000e+00i
      ];
    
      abcd2 = [
        0.000000000000000e+00 + 0.000000000000000e+00i
        9.996875162757026e-01 + 2.499739591471233e-02i
        0.000000000000000e+00 + 0.000000000000000e+00i
        -9.637765873613319e-01 - 2.667108727707685e-01i
      ];
  
      abcd3 = [
        9.987502603949663e-01 - 4.997916927067833e-02i
        1.000000000000000e+00 + 0.000000000000000e+00i
        0.000000000000000e+00 + 0.000000000000000e+00i
        0.000000000000000e+00 + 0.000000000000000e+00i
      ];
    
      if hat
        q13 = [1,2];
        q2 = [0,1];
      else
        q13 = [0,1];
        q2 = [1,2];
      end
      
      gate1 = f3c.qgates.RotationTFXYMatrix( q13, abcd1 );
      gate2 = f3c.qgates.RotationTFXYMatrix( q2, abcd2 );
      gate3 = f3c.qgates.RotationTFXYMatrix( q13, abcd3 );
      test.runTurnoverTest( gate1, gate2, gate3, tol, test ) ;
      % ------------------------------------------------------------------------
     
      hat = true;
      abcd1 = [
        3.987390988927981e-01 - 9.170644094141694e-01i
        -8.851109971106588e-01 + 4.653799767864700e-01i
        3.767165746647338e-11 - 1.436744797632233e-10i
        -1.477308631374141e-10 - 1.539298066265034e-11i
      ];
    
      abcd2 = [
        9.754613486067534e-03 - 7.745308712318724e-03i
        1.378034074544130e-02 + 0.000000000000000e+00i
        -2.159081679304937e-01 - 9.763164267950417e-01i
        -9.999224258454700e-01 + 8.131516293641283e-20i
      ];
  
      abcd3 = [
        -1.049695306579096e-02 + 6.060574579217108e-03i
        1.212091530766865e-02 + 1.292469707114106e-26i
        -7.795229670875499e-01 + 6.262563590057401e-01i
        -9.999265390077934e-01 - 1.387778780781446e-17i
      ];
    
      if hat
        q13 = [1,2];
        q2 = [0,1];
      else
        q13 = [0,1];
        q2 = [1,2];
      end
      
      gate1 = f3c.qgates.RotationTFXYMatrix( q13, abcd1 );
      gate2 = f3c.qgates.RotationTFXYMatrix( q2, abcd2 );
      gate3 = f3c.qgates.RotationTFXYMatrix( q13, abcd3 );
      test.runTurnoverTest( gate1, gate2, gate3, tol, test ) ;
      % ------------------------------------------------------------------------
    end
  end
  
  methods (Static)
    
    function runTurnoverTest( gate1, gate2, gate3, tol, test )
       % turnover
       [gateA, gateB, gateC] = turnover( gate1, gate2, gate3 );
       
       % Cin
       Cin = qclab.QCircuit( 3 );
       Cin.push_back( gate1 );
       Cin.push_back( gate2 );
       Cin.push_back( gate3 );
       matin = Cin.matrix;
       
       
       % Cout
       Cout = qclab.QCircuit( 3 ) ;
       Cout.push_back( gateA );
       Cout.push_back( gateB );
       Cout.push_back( gateC );
       
       matout = Cout.matrix;
       
       error = norm( matin - matout );
       
       if error > tol
         fprintf(1,'error:\t%e\n', error );
         fprintf(1,'gate 1:\n');
         v = gate1.abcd;
         v(:)
         fprintf(1,'gate 2:\n');
         v = gate2.abcd;
         v(:)
         fprintf(1,'gate 3:\n');
         v = gate3.abcd;
         v(:)
       end
       
       test.verifyEqual( error, 0, 'AbsTol', tol );
    end
    
    function runEdgeTurnoverTest( ratios, hat, tol, test )
      if hat
        gate1 = f3c.qgates.RotationTFXYMatrix( [1,2], [1 1 0 0] );
        gate2 = f3c.qgates.RotationTFXYMatrix( [0,1], [1 1 0 0] );
        gate3 = f3c.qgates.RotationTFXYMatrix( [1,2], [1 1 0 0] );
      else
        gate1 = f3c.qgates.RotationTFXYMatrix( [0,1], [1 1 0 0] );
        gate2 = f3c.qgates.RotationTFXYMatrix( [1,2], [1 1 0 0] );
        gate3 = f3c.qgates.RotationTFXYMatrix( [0,1], [1 1 0 0] );
      end
      
      for i1 = 1:2
        for i2 = 1:2
          for i3 = 1:2
            for i4 = 1:2
              for i5 = 1:2
                for i6 = 1:2
                  
                  [a1, d1] = test.genSU2mat(ratios(i1));
                  [b1, c1] = test.genSU2mat(ratios(i2));
                  gate1.update([a1, b1, c1, d1]);
                  
                  
                  [a2, d2] = test.genSU2mat(ratios(i3));
                  [b2, c2] = test.genSU2mat(ratios(i4));
                  gate2.update([a2, b2, c2, d2]);
                  
                  [a3, d3] = test.genSU2mat(ratios(i5));
                  [b3, c3] = test.genSU2mat(ratios(i6));
                  gate3.update([a3, b3, c3, d3]);
                  
                  test.runTurnoverTest( gate1, gate2, gate3, tol, test )
                 
                  
                end
              end
            end
          end
        end
      end      
    end
    
    function [alpha, beta] = genSU2mat(DtoAD)
      beta = 1 / sqrt(1 + DtoAD^2);
      alpha = sqrt(1 - beta^2);
      alpha = alpha * exp( 1i * randn );
      beta = beta * exp( 1i * randn );
    end
  end
end