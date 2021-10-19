%% OVERVIEW
% ======================================================================= %
%{ 
    This program simulates time crystalline order emergent from a Floquet 
    disordered system of spins, by observing the  dynamics 2-time 
    correlator.
    
    REFERENCES:
        - Xiao Mi, et al. arxiv:2107.13571, 2021.
        - Matteo Ippoliti, et al. PRX Quantum 2(3):030346, 2021.
%}
% ======================================================================= %

% ======================================================================= %
%{ 
    author: REUBEN R. W. WANG 
    email: reuben.wang@colorado.edu 
%}
% ======================================================================= %

warning('off', 'all')
%% INITIALIZING THE PHYSICAL SYSTEM
% ======================================================================= %
% system parameters
N = 10; % number of spin sites
J0 = 1; % strength of interaction
g = 0.94; % pi-pulse distortion parameter
nT = 100; % number of periodic drives to evolve by
% ======================================================================= %

% constructing the Floquet unitary operators
% ======================================================================= %
UF = FloquetUnitary (N, g, J0);
UF_dagger = UF';

% diagonalizing Floquet unitary
[V,D] = eigs(UF, 2^N);
V = sparse(V); D = sparse(D);
disp('Floquet unitary has been diagonalized...')
% ======================================================================= %

%% RUNNING DYNAMICAL EVOLUTION
% ======================================================================= %
%{
    OPTIONAL ARGUMENTS
        - 'data2txt': set to false to turn off text writing.
        - 'spinSite': specify site of the observable spin.
        - 'initState': choose init. state from {'Neel', 'random'}.
%}
output = timeCrystal (N, g, nT, V, D, 'spinSite', 3);
time = 1:(nT+1);
expectZt = output{1}; 
expectZ0Zt = output{2}; 
% ======================================================================= %

%% PLOTTING RESULTS
% ======================================================================= %
% plotting <Z(t)>
fig1 = figure;
axes1 = axes('Parent', fig1);
hold(axes1,'on');
plot(time, expectZt, 'b--', 'LineWidth', 1.0)
scatter(time, expectZt, 25, 'b', 'filled')
hold(axes1,'off');
xlabel('$t/T$', 'Interpreter', 'latex')
ylabel('$\langle Z(t) \rangle$', 'Interpreter', 'latex')
title(['no. of sites = ', num2str(N), ', $g =$ ', num2str(g),...
       ', $J_0 = $ ', num2str(J0),'$\pi$'], 'Interpreter', 'latex')
   
% plotting <Z(0)Z(t)>
fig2 = figure;
axes2 = axes('Parent', fig2);
hold(axes2,'on');
scatter(time, zeros(size(time)), 25, 'r', 'filled')
plot(time, expectZ0Zt, 'b--', 'LineWidth', 1.0)
scatter(time, expectZ0Zt, 25, 'b', 'filled')
xlabel('$t/T$', 'Interpreter', 'latex')
ylabel('$\langle Z(0) Z(t) \rangle$', 'Interpreter', 'latex')
title(['no. of sites = ', num2str(N), ', $g =$ ', num2str(g),...
       ', $J_0 = $ ', num2str(J0),'$\pi$'], 'Interpreter', 'latex')
% ======================================================================= %
%%




% ======================================================================= %
% ======================================================================= %
% ======================================================================= %
% ============================= BLANK SPACE ============================= %
% ======================================================================= %
% ======================================================================= %
% ======================================================================= %





%% FUNCTION DEFINITIONS
% ======================================================================= %
function output = timeCrystal (N, g, nT, V, D, varargin)
% ======================================================================= %
%{
    Description
        Method that constructs a Pauli at site site_idx in the full 
        Hilbert space representation.

    Arguments
        - N: no. of spin sites.
        - g: pi-pulse parameter.
        - nT: number of periods to evolve system.
        - V, D: eigenvectors and eigenvalues of UF.
        - varargin: optional arguments.
%}
% ======================================================================= %

% allow optional inputs of arguments
% ======================================================================= %
% addParamValue
IP = inputParser; 

% optional arguments
IP.addParameter('data2txt',true); % choose to write data to txt file
IP.addParameter('spinSite',2); % site of the observable spin
IP.addParameter('initState','Neel'); % choose init. state of system

% parsing arguments
IP.parse(varargin{:});
% ======================================================================= %

% error checks
% ======================================================================= %
if (N <= 2)
error('ERROR: Too few spins are being simulated! Must set N > 2.') 
end
% ======================================================================= %

% initial state preparation
% ======================================================================= %
obsSite = IP.Results.spinSite;
if strcmp(IP.Results.initState, 'Neel')
state0 = NeelState (N); % init. Neel state
elseif strcmp(IP.Results.initState, 'random')
state0 = round(rand(2^N, 1)); % init. random state 
end
state0 = sparse(state0/vecnorm(state0)); % normalize state
% ======================================================================= %

% init. observable
% ======================================================================= %
% identity and Pauli matrices
Id = sparse(eye(2));
% sigX = sparse([0, 1; 1, 0]);
% sigY = sparse([0, -1j; 1j, 0]);
sigZ = sparse([1, 0; 0 -1]);

% obs. spin operators
obsZ = sparse(oneSitePauli(N, obsSite, Id, Id, sigZ)); 
% obsSite_2 = 7; % site of observed spin
% obsZ_2 = sparse(oneSitePauli(N, obsSite_2, Id, Id, sigZ)); % obs. spin operator
% ======================================================================= %

% computing <Z(t)>, <Z(0)Z(t)> and <Z_i(t)Z_j(t)>
% ======================================================================= %
expectZt = zeros(1, nT+1); % init. <Z(t)> array
expectZt(1) = state0' * obsZ * state0; % init. <Z(0)> 
expectZ0Zt = zeros(1, nT+1); % init. <Z(0)Z(t)> array
expectZ0Zt(1) = state0' * obsZ^2 * state0; % init. <Z(0)Z(0)> 
% expectZiZj = zeros(1, nT+1); % init. <Z_i(t)Z_j(t)> array
% expectZiZj(1) = state0' * (obsZ*obsZ_2) * state0; % init. <Z(0)Z(0)>
tic; % init. timer
for tstep = 2:(nT+1)  
    % unitary matrix at current time
    UFt = sparse(V*(D.^(tstep-1))/V); 
    state = sparse(UFt * state0);
    
    % expectation values
    expectZt(tstep) = state' * obsZ * state;
    expectZ0Zt(tstep) = state0' * obsZ * (UFt' * obsZ) * state;
    % expectZiZj(tstep) = state' * (obsZ*obsZ_2) * state;
    
    % close timer
    if (tstep == 6)
    end_time = toc;
    disp(['etimated runtime: ', num2str(end_time*nT/5),' (s)'])
    end
end
toc % print actual runtime
% ======================================================================= %

% storing data as text file
% ======================================================================= %
if (IP.Results.data2txt)
time = 1:(nT+1);
write2txt (time, expectZt, expectZ0Zt, N, g);
end
% ======================================================================= %

% return function outputs
% ======================================================================= %
output = {expectZt, expectZ0Zt};
% ======================================================================= %
end

function write2txt (time, expectZt, expectZ0Zt, N, g, varargin)
% ======================================================================= %
%{
    Description
        Method that writes the input data to a text file.

    Arguments
        - time: time array. 
        - expectZt: <Z(t)> array.
        - expectZ0Zt: <Z(0)Z(t)> array.
        - N, g: no. spins and pi-pulse parameter.
        - varargin: optional arguments.
%}
% ======================================================================= %

% allow optional inputs of arguments
% ======================================================================= %
% addParamValue
IP = inputParser; 

% optional arguments
IP.addParameter('dataPath',pwd); % data path

% parsing arguments
IP.parse(varargin{:});
% ======================================================================= %

% writing data to text file
% ======================================================================= %
writeMat = ["t/T", "<Z(t)>", "<Z(0)Z(t)>"; 
            join(string(num2str(time(:), '%0.5e')), 2),...
            join(string(num2str(expectZt(:), '%0.5e')), 2),...
            join(string(num2str(expectZ0Zt(:), '%0.5e')), 2)];
fileName = strcat(strcat(strcat("N=", num2str(N)),...
                         strcat("_g=", num2str(g*100))),...
                  "e-2.txt");
curretFolder = IP.Results.dataPath;
curretFolder = strcat(curretFolder, '/');
filePath = strcat(curretFolder, fileName);
disp(' '); disp('text file being written to:'); disp(filePath)
writematrix(writeMat, filePath, 'Delimiter','tab');
% ======================================================================= %
end

function UF = FloquetUnitary (N, g, J0)
% ======================================================================= %
%{
    Description
        Method that computes the Floquet unitary matrix.

    Arguments
        - N: no. of spins.
        - g: pi-pulse parameter.
        - J0: spin-spin interaction strength.
%}
% ======================================================================= %

% identity and Pauli matrices
% ======================================================================= %
Id = sparse(eye(2));
sigX = sparse([0, 1; 1, 0]);
% sigY = sparse([0, -1j; 1j, 0]);
sigZ = sparse([1, 0; 0 -1]);
% ======================================================================= %

% logitudinal field unitary
% ======================================================================= %
b_array = -pi + 2*pi*rand(1, N); % random parameters
H_long = 0; % logituddinal field Hamiltonian
for site_idx = 1:N
    if (site_idx == 1)
        H_long = H_long +...
            b_array(site_idx)*oneSitePauli(N, site_idx, sigZ, Id, sigZ); 
    else
        H_long = H_long +...
            b_array(site_idx)*oneSitePauli(N, site_idx, Id, Id, sigZ); 
    end
end
UF_long = sparse(expm(-1j/2 * H_long)); 
% ======================================================================= %

% Ising-interaction unitary
% ======================================================================= %
J_array = -J0*pi + (-0.5 + rand(1, N-1))*pi; % random parameters
H_NNZ = 0; % X-rotation Hamiltonian
for site_idx = 1:(N-1)
    if (site_idx == 1)
        H_NNZ = H_NNZ +...
            J_array(site_idx)*NNPauli(N, site_idx, sigZ, Id, sigZ); 
    else
        H_NNZ = H_NNZ +...
            J_array(site_idx)*NNPauli(N, site_idx, Id, Id, sigZ); 
    end
end
UF_NNZ = sparse(expm(-1j/4 * H_NNZ)); 
% ======================================================================= %

% X-rotation unitary
% ======================================================================= %
H_Xrot = 0; % X-rotation Hamiltonian
for site_idx = 1:N
    if (site_idx == 1)
        H_Xrot = H_Xrot + oneSitePauli(N, site_idx, sigX, Id, sigX); 
    else
        H_Xrot = H_Xrot + oneSitePauli(N, site_idx, Id, Id, sigX); 
    end
end
UF_Xrot = sparse(expm(-1j/2 * pi*g * H_Xrot)); 
% ======================================================================= %

% total Floquet unitary
% ======================================================================= %
UF = UF_long * UF_NNZ * UF_Xrot;
% ======================================================================= %
end

function outMat = oneSitePauli (N, site_idx, M, Id, Pauli)
% ======================================================================= %
%{
    Description
        Method that constructs a Pauli at site site_idx in the full 
        Hilbert space representation.

    Arguments
        - N: no. of spin sites.
        - site_idx: spin site index.
        - M: matrix of current iteration.
        - I, Pauli: identity and Pauli matrices.
%}
% ======================================================================= %

% running recursion iteration
% ======================================================================= %
if (N <= 1)
    outMat = M;
else
    if (log2(length(M))+1 == site_idx)
        outMat = oneSitePauli (N-1, site_idx, kron(M, Pauli), Id, Pauli);
    else
        outMat = oneSitePauli (N-1, site_idx, kron(M, Id), Id, Pauli);
    end
end
% ======================================================================= %
end

function outMat = NNPauli (N, site_idx, M, Id, Pauli)
% ======================================================================= %
%{
    Description
        Method that constructs a nearest neighbor Pauli interactions at 
        site site_idx and site_idx+1 in the full Hilbert space 
        representation.

    Arguments
        - N: no. of spin sites.
        - site_idx: spin site index.
        - M: matrix of current iteration.
        - I, Pauli: identity and Pauli matrices.
%}
% ======================================================================= %

% running recursion iteration
% ======================================================================= %
if (N <= 1)
    outMat = M;
else
    if (log2(length(M))+1 == site_idx) || (log2(length(M))+1 == site_idx+1)
        outMat = NNPauli (N-1, site_idx, kron(M, Pauli), Id, Pauli);
    else
        outMat = NNPauli (N-1, site_idx, kron(M, Id), Id, Pauli);
    end
end
% ======================================================================= %
end

function outState = NeelState (N)
% ======================================================================= %
%{
    Description
        Method that constructs a Neel state.

    Arguments
        - N: no. of spin sites.
%}
% ======================================================================= %

% constructing the Neel state
% ======================================================================= %
% defining one-site qubit states
zero = [1; 0];
one = [0; 1];
pair = kron(zero, one);

% taking tensor products
if (mod(N,2) == 0) % even no. of sites
    outState = recursiveKron ( N/2, pair, pair );
else
    nm1_state = recursiveKron ( (N-1)/2, pair, pair );
    outState = kron(nm1_state, zero);
end
% ======================================================================= %
end

function outMat = recursiveKron ( num, A0, A )
% ======================================================================= %
%{
    Description
        Method that constructs a nearest neighbor Pauli interactions at 
        site site_idx and site_idx+1 in the full Hilbert space 
        representation.

    Arguments
        - num: no. of times to run recursion.
        - A0: matrix to be tensored recursively.
        - A: current matrix in recursion.
%}
% ======================================================================= %

% running recursion iteration
% ======================================================================= %
if (num <= 1)
    outMat = A;
else
    outMat = recursiveKron ( num-1, A0, kron(A0, A) );
end
% ======================================================================= %
end
% ======================================================================= %
