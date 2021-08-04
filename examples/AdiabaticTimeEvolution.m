% ------------------------------------------------------------------------------
% Adiabatic Time Evolution for the TFXY Hamiltonian. The parameters of this
% script are setup to compute the circuits used for the adiabatic state
% preparation experiment presented in Fig. 6 of:
%     Algebraic Compression of Quantum Circuits for Hamiltonian Evolution, 
%     Efekan Kockcu, Daan Camps, Lindsay Bassman, J. K. Freericks, 
%     Wibe A. de Jong, Roel Van Beeumen, and Alexander F. Kemper (2021) 
% ------------------------------------------------------------------------------

% ### Parameters to adjust #####################################################

% Overall parameters -----------------------------------------------------------
N = 5;                             % number of spins in the TFXY Hamiltonian
nsteps = 240;                      % number of Trotter steps: 
                                   %  * 240 for simulation with large Trotter
                                   %    error (purple line in Fig. 6)
                                   %  * 1200 for simulation with small Trotter
                                   %    error (red line in Fig. 6)
dt = 0.25;                         % size of a Trotter step:
                                   %  * 0.25 for simulation with large Trotter
                                   %    error (purple line in Fig. 6)
                                   %  * 0.05 for simulation with small Trotter
                                   %    error (purple line in Fig. 6)
% circuit storage directory
folder = 'ADD LOCAL PATH HERE';
saveCircuit = false;               % true/false: save quantum circuits to file.
% ------------------------------------------------------------------------------

% Transverse field hz ----------------------------------------------------------
hz.start = -1.0;                      % value of hz at start of simulation
hz.end = -1.0;                        % value of hz at end of simulation
hz.startRamp = 1;                     % step when ramping starts
hz.endRamp = nsteps;                  % step when ramping ends
hz.rampType = @linearRamp;            % type of ramp: linear or sigmoid
% ------------------------------------------------------------------------------

% Coupling parameter Jx --------------------------------------------------------
Jx.start = 0.0;                       % value of Jx at start of simulation
Jx.end = -2.0;                        % value of Jx at end of simulation
Jx.startRamp = 1;                     % step when ramping starts
Jx.endRamp = round( 0.5 * nsteps );   % step when ramping ends
Jx.rampType = @linearRamp;            % type of ramp: linear or sigmoid
% ------------------------------------------------------------------------------

% Coupling parameter Jy --------------------------------------------------------
Jy.start = 0.0;                       % value of Jy at start of simulation
Jy.end = 0.0;                         % value of Jy at end of simulation
Jy.startRamp = 1;                     % step when ramping starts
Jy.endRamp = nsteps;                  % step when ramping ends
Jy.rampType = @linearRamp;            % type of ramp: linear or sigmoid
% ------------------------------------------------------------------------------

% ############################################################################## 
clc
fprintf( 1, '\n---------------------------------------------------------------\n');
fprintf( 1, 'Generating and validating quantum circuits for TFXY Hamiltonian\n' );
fprintf( 1, '---------------------------------------------------------------\n\n');

% generate all parameter arrays for full simulation.
hz_arr = arrayFromCell( hz, nsteps );
Jx_arr = arrayFromCell( Jx, nsteps );
Jy_arr = arrayFromCell( Jy, nsteps );

% set timing arrays.
tfinal = ( nsteps - 1 ) * dt;
t_arr = linspace( 0, tfinal, nsteps );

% print information about simulation to command window.
fprintf( 1, 'Number of spins = \t\t\t%d\n', N );
fprintf( 1, 'Number of Trotter steps = \t\t%d\n', nsteps );
fprintf( 1, 'Stepsize of Trotter step = \t\t%.3f\n', dt );
fprintf( 1, 'Final simulation time = \t\t%.3f\n', tfinal );

% save some plot colors.
color1 = [252, 186, 3]./255;
color2 = [3, 92, 59]./255;
color3 = [167, 38, 222]./255;
color4 = [37, 122, 118]./255;

% plot the simulation parameter curves.
figure(1); clf
plot(t_arr, Jx_arr, 'color', color1, 'LineWidth', 2, 'DisplayName', 'J_x');
hold on
plot(t_arr, Jy_arr, 'color', color2, 'LineWidth', 2,'DisplayName', 'J_y');
plot(t_arr, hz_arr, 'color', color3, 'LineWidth', 2,'DisplayName', 'h_z');
xlabel('Simulation time');
title('Simulation parameters for TFXY');
legend();

% compute the Hamiltonian eigenvalues throughout simulation.
energies = zeros(nsteps, 2^N);
for i = 1:nsteps
  energies(i, :) = eig( hamTFXY(N, Jx_arr(i), Jy_arr(i), hz_arr(i)) );
end

% compute the magnetization observable O
Z = qclab.qgates.PauliZ.matrix;
O = zeros(2^N, 2^N);
for  i = 0:N-1
  O = O + kron(kron( qclab.qId(i), Z), qclab.qId(N - i - 1) );
end

% start state(s) for simulation.
psi0 = zeros(2^N, 1);
psi0(end) = 1; % all 1's initial state.
psi_ex = psi0;
psi_tr = psi0;

% setup arrays for results.

% energy
energy_ex = zeros( nsteps, 1 );  % energy with exact unitary
energy_tr = zeros( nsteps, 1 );  % energy with trotter unitary
energy_fc = zeros( nsteps, 1 );  % energy with full circuit
energy_cc = zeros( nsteps, 1 );  % energy with compressed circuit

% magnetization
mag_ex = zeros( nsteps, 1 ); % magnetization with exact unitary
mag_tr = zeros( nsteps, 1 ); % magnetization with trotter unitary
mag_fc = zeros( nsteps, 1 ); % magnetization with full circuit
mag_cc = zeros( nsteps, 1 ); % magnetization with compressed circuit

% do the full simulation and generate + save circuits.
fullCircuitMat = qclab.QCircuit( N );
fullCircuitRot = qclab.QCircuit( N );
compressedCircuitRot = qclab.QCircuit( N );
TFXY = @f3c.qgates.RotationTFXY;

numerical_error = zeros( nsteps, 1 ); % numerical error on compression algorithm

circLayers = 1;
for i = 1:nsteps
  % Print info to command window.
  fprintf( 1, 'Current time = \t\t\t%.3f\t(%d/%d)\n', t_arr(i), i, nsteps );
  % Hamiltonian at this timestep.
  H = hamTFXY( N, Jx_arr(i), Jy_arr(i), hz_arr(i) );
  
  % full unitaries: exact and trotter ------------------------------------------
  psi_ex = expm( -1i * dt *  H ) * psi_ex ;
  psi_tr = expm( -1i * dt * hamTFXY( N, Jx_arr(i), Jy_arr(i), 0 ) ) * ...
           expm( -1i * dt * hamTFXY( N, 0, 0, hz_arr(i) ) ) * psi_tr ;
  % energies
  energy_ex(i) = real( psi_ex' * H * psi_ex );
  energy_tr(i) = real( psi_tr' * H * psi_tr );
  % magnetization
  mag_ex(i) = -real( psi_ex' * O * psi_ex )/N;
  mag_tr(i) = -real( psi_tr' * O * psi_tr )/N;
  
  % quantum circuits -----------------------------------------------------------
  [timeStepMat, timeStepRot] = f3c.timeStepTFXY( N, dt, Jx_arr(i), Jy_arr(i), hz_arr(i) );
  if circLayers < ceil(N/2) % not yet at minimal depth size
    for j = 1:timeStepMat.nbGates
      fullCircuitMat.push_back( timeStepMat.gates(j) );
      fullCircuitRot.push_back( timeStepRot.gates(j) );
    end
    squareCircuit = copy( fullCircuitMat );
  elseif circLayers == ceil(N/2) % becomes minimal depth at current size
    if mod( N, 2 ) == 0 % even number of spins
      for j = 1:timeStepMat.nbGates
        fullCircuitMat.push_back( timeStepMat.gates(j) );
        fullCircuitRot.push_back( timeStepRot.gates(j) );
      end
      squareCircuit = f3c.SquareCircuit( fullCircuitMat );
      triangleCircuit = f3c.TriangleCircuit( squareCircuit );
    else
      % add first vertical layer to make it square
      for j = 1:(N-1)/2
        fullCircuitMat.push_back( timeStepMat.gates(j) );
        fullCircuitRot.push_back( timeStepRot.gates(j) );
      end
      squareCircuit = f3c.SquareCircuit( fullCircuitMat );
      triangleCircuit = f3c.TriangleCircuit( squareCircuit );
      % add second vertical layer to fullCircuit and merge with triangle
      for j = ((N-1)/2)+1:timeStepMat.nbGates
        fullCircuitMat.push_back( timeStepMat.gates(j) );
        triangleCircuit.merge( timeStepMat.gates(j), 'R' );
        fullCircuitRot.push_back( timeStepRot.gates(j) );
      end
      squareCircuit = f3c.SquareCircuit( triangleCircuit );
    end
  else % larger than minimal depth
    for j = 1:timeStepMat.nbGates
      fullCircuitMat.push_back( timeStepMat.gates(j) );
      triangleCircuit.merge( timeStepMat.gates(j), 'R' );
      
      fullCircuitRot.push_back( timeStepRot.gates(j) );
    end
    squareCircuit = f3c.SquareCircuit( triangleCircuit );
  end
  
  % rotation representation of compressed circuit.
  compressedCircuitRot.clear();
  for j = 1:squareCircuit.nbGates
    compressedCircuitRot.push_back( TFXY( squareCircuit.gates(j) ) );
  end
  
  % compute the state based on the quantum circuits
  psi_fc = fullCircuitMat.apply('R', 'N', N, psi0 );
  psi_cc = squareCircuit.apply('R', 'N', N, psi0 );
  
  % energies
  energy_fc(i) = real( psi_fc' * H * psi_fc );
  energy_cc(i) = real( psi_cc' * H * psi_cc );
  % magnetization
  mag_fc(i) = -real( psi_fc' * O * psi_fc )/N;
  mag_cc(i) = -real( psi_cc' * O * psi_cc )/N;
  
  numerical_error(i) = ...
    norm( fullCircuitRot.matrix - compressedCircuitRot.matrix, 'fro' );
  if mod(i, 10) == 0
    fprintf(1, 'Numerical error on compressed circuit: %e\n', ...
                                                           numerical_error(i) );
  end
  % save the circuits to file.
  if saveCircuit
    fname_fc = [folder,'TFXY_full_circuit_t_',num2str(t_arr(i)),'.qasm'];
    fid = fopen( fname_fc, 'w+' );
    fprintf( fid, '// Generated by f3c\n');
    fprintf( fid, '// https://github.com/QuantumComputingLab/f3c\n\n' );
    fprintf( fid, '// Current time = \t%6.4f [dt = %6.4f][timestep = %d/%d]\n', ...
                                                      t_arr(i), dt, i, nsteps );
    fprintf( fid, '// Jx = %6.4f,\t Jy = %6.4f,\t hz = %6.4f\n\n', ...
                                              Jx_arr(i), Jy_arr(i), hz_arr(i) );
    fprintf( fid, 'OPENQASM 2.0;\ninclude "qelib1.inc";\n\n');
    fprintf( fid, 'qreg q[%d];\n',N);
    fullCircuitRot.toQASM( fid );
    fclose( fid );
    
    fname_cc = [folder,'TFXY_compressed_circuit_t_',num2str(t_arr(i)),'.qasm'];
    fid = fopen( fname_cc, 'w+' );
    fprintf( fid, '// Generated by f3c\n');
    fprintf( fid, '// https://github.com/QuantumComputingLab/f3c\n\n' );
    fprintf( fid, '// Current time = \t%6.4f [dt = %6.4f][timestep = %d/%d]\n', ...
                                                      t_arr(i), dt, i, nsteps );
    fprintf( fid, '// Jx = %6.4f,\t Jy = %6.4f,\t hz = %6.4f\n\n', ...
                                              Jx_arr(i), Jy_arr(i), hz_arr(i) );
    fprintf( fid, 'OPENQASM 2.0;\ninclude "qelib1.inc";\n\n');
    fprintf( fid, 'qreg q[%d];\n',N);
    compressedCircuitRot.toQASM( fid );
    fclose( fid );
  end 
  
  circLayers = circLayers + 1;
end

% plot the Hamiltonian eigenvalues throughout simulation together
% with approximate ground states.
figure(2); clf
axEnergies = axes;
p = plot(t_arr, energies, 'k-');
for i = 1:length(p)
    p(i).Annotation.LegendInformation.IconDisplayStyle = 'off';
end
hold on
plot( t_arr, energy_ex, 'color', color1, 'LineWidth', 2, 'DisplayName', 'Exact Unitary' );
plot( t_arr, energy_tr, 'color', color2, 'LineWidth', 2, 'DisplayName', 'Trotter Unitary');
plot( t_arr, energy_fc, '--','color', color3, 'LineWidth', 2, 'DisplayName', 'Full Circuit');
plot( t_arr, energy_cc, ':','color', color4, 'LineWidth', 2, 'DisplayName', 'Compressed Circuit');
xlabel('Simulation time');
title('Energies during TFXY simulation');
legend();

% plot the magnetization.
figure(3); clf
plot( t_arr, mag_ex, 'color', color1, 'LineWidth', 2, 'DisplayName', 'Exact Unitary' );
hold on
plot( t_arr, mag_tr, 'color', color2, 'LineWidth', 2, 'DisplayName', 'Trotter Unitary');
plot( t_arr, mag_fc, '--','color', color3, 'LineWidth', 2, 'DisplayName', 'Full Circuit');
plot( t_arr, mag_cc, ':','color', color4, 'LineWidth', 2, 'DisplayName', 'Compressed Circuit');
xlabel('Simulation time');
title('Magnetization during TFXY simulation');
legend();


% plot the numerical error on the compression algorithm
figure(4); clf
loglog( t_arr, numerical_error, 'color', color1, 'LineWidth', 2 );
xlabel('Simulation time');
title('Numerical error on compression algorithm');
ylabel('Frobenius norm of error');
ylim([1e-16,1e0])

% -- Functions ---------------------------------------------------------------
function f = arrayFromCell( in, nsteps )
  f = in.rampType( in.startRamp, in.endRamp, nsteps );
  f = in.start + f .* ( in.end - in.start );
end

% generates a linear ramp from 0 to 1 with npoints.
% it starts ramping up from point start \in [1, npoints-1]
% and ends ramping up at point end \in [2, npoints]
function f = linearRamp( startpoint, endpoint, npoints )
  f = zeros( npoints, 1 );
  ramp = linspace( 0, 1, endpoint - startpoint + 1);
  f(startpoint:endpoint) = ramp;
  f(endpoint+1:end) = 1;
end

% generates a sigmoid ramp from 0 to 1 with npoints
% it starts ramping up from point start \in [1, npoints-1]
% and ends ramping up at point end \in [2, npoints]
function f = sigmoidRamp( startpoint, endpoint, npoints )
  f = zeros( npoints, 1 );
  x = linspace( -10, 10, endpoint - startpoint + 1 );
  f(startpoint:endpoint) = 1./( 1 + exp(-x) );
  f(endpoint+1:end) = 1;  
end

% generate the TFXY Hamiltonian for ordered (scalar) and disordered models
% (arrays)
function H = hamTFXY(N, Jx, Jy, hz)
  assert( qclab.isNonNegInteger( N ) );
  if isscalar( Jx )
    Jx = repmat( Jx, N-1, 1 );
  end
  assert( length(Jx) == N-1 );
  if isscalar( Jy )
    Jy = repmat( Jy, N-1, 1 );
  end
  assert( length(Jy) == N-1 );
  if isscalar( hz )
    hz = repmat( hz, N, 1 );
  end
  assert( length(hz) == N );
  XX = kron( qclab.qgates.PauliX.matrix, qclab.qgates.PauliX.matrix );
  YY = kron( qclab.qgates.PauliY.matrix, qclab.qgates.PauliY.matrix );
  Z = qclab.qgates.PauliZ.matrix;
  H = zeros(2^N,2^N);
  for i=0:N-2
    H = H + Jx(i+1) * kron(kron( qclab.qId(i), XX ), qclab.qId(N - i - 2) ) ...
          + Jy(i+1) * kron(kron( qclab.qId(i), YY ), qclab.qId(N - i - 2) );
  end
  for i=0:N-1
    H = H + hz(i+1) * kron(kron( qclab.qId(i), Z), qclab.qId(N - i - 1) );
  end
end