 %clc;close all; 
load('R-PAM_DL_DO_150_M8.mat','C')    % Load the codebook

M=size(C,2);  % size of the codebook for one user
m=log2(M); % no of bits in a symbol/block
K=4;   % No of orthogonal PREs/ OFDM subcarriers
J=6;   % No of users/layers
Nr=8;Nt=1;
F=get_indicator_matrix_factor_graph(C, J, K);
F1=zeros(Nr*K,Nt*J);
for nr=1:Nr
    for nt=1:Nt
    
F1((nr-1)*K+1:nr*K,(nt-1)*J+1:nt*J)=F;
    end
    
end
%% %%%%%%% Power Allocation %%%%%%%%%
% power_users=ones(1,J);  
% sqrt_power_matrix=ones(K,J);
%% Monte Carlo %%%%%%%%%%%%%%%%%%%%%%%%
%Eb_N0_dB=5:3:30;
Eb_N0_dB=-5:2:5;
%Eb_N0_dB=-5:2:5;
Eb_N0=10.^(Eb_N0_dB/10);
Es=sum(sum((abs(C)).^2))/length(C);
Eb=Es/m;   
N0=Eb./Eb_N0;
sigma=sqrt(N0/2);
max_block_errors=[   400 350 250 200 150 100];
%max_block_errors= [2 1 1]  ;
%max_block_errors=1;
T=10;  % maximum number of iterations for SCMA detection
SER=zeros(1,length(Eb_N0_dB));
%cf=500; % LLR clipping factor (used only in the LLR-domain MPA)
%% 
 for e=1:length(Eb_N0_dB)

    block_errors=0;
    symbol_errors=0;
    block=0;
    while block_errors<max_block_errors(e)
        %%   SCMA Encoding %%%%%%
        bits_A1=randi([0 1],J,m);% blocks of bits for all users Tx-ant-1
%         bits_A2=randi([0 1],J,m);% blocks of bits for all users Tx-ant-2
%         bits_A3=randi([0 1],J,m);% blocks of bits for all users Tx-ant-3
%         bits_A4=randi([0 1],J,m);% blocks of bits for all users Tx-ant-4
%         bits_A5=randi([0 1],J,m);% blocks of bits for all users Tx-ant-1
%         bits_A6=randi([0 1],J,m);% blocks of bits for all users Tx-ant-2
%         bits_A7=randi([0 1],J,m);% blocks of bits for all users Tx-ant-3
%         bits_A8=randi([0 1],J,m);% blocks of bits for all users Tx-ant-4
        %mapping
        symbol_indices_A1=bi2de(bits_A1,'left-msb')+1; % symbols for all users Tx-1
%         symbol_indices_A2=bi2de(bits_A2,'left-msb')+1; % symbols for all users Tx-2
%         symbol_indices_A3=bi2de(bits_A3,'left-msb')+1; % symbols for all users Tx-3
%         symbol_indices_A4=bi2de(bits_A4,'left-msb')+1; % symbols for all users Tx-2
%         
%         symbol_indices_A5=bi2de(bits_A5,'left-msb')+1; % symbols for all users Tx-1
%         symbol_indices_A6=bi2de(bits_A6,'left-msb')+1; % symbols for all users Tx-2
%         symbol_indices_A7=bi2de(bits_A7,'left-msb')+1; % symbols for all users Tx-3
%         symbol_indices_A8=bi2de(bits_A8,'left-msb')+1; % symbols for all users Tx-2
        %SCMA_codewords=zeros(K,J);  % collection of the codewords for all users
        SCMA_codewords_A1=zeros(K,J);  % collection of the codewords for all users Tx-1 
%        SCMA_codewords_A2=zeros(K,J);  % collection of the codewords for all users Tx-2
%         SCMA_codewords_A3=zeros(K,J);  % collection of the codewords for all users Tx-3
%         SCMA_codewords_A4=zeros(K,J);  % collection of the codewords for all users Tx-4 
%         
%         SCMA_codewords_A5=zeros(K,J);  % collection of the codewords for all users Tx-1 
%         SCMA_codewords_A6=zeros(K,J);  % collection of the codewords for all users Tx-2
%         SCMA_codewords_A7=zeros(K,J);  % collection of the codewords for all users Tx-3
%         SCMA_codewords_A8=zeros(K,J); 
        for j=1:J         % for each user
            present_codebook=C((j-1)*K+1:j*K,:);   % codebook for the jth user
            SCMA_codewords_A1(:,j)=present_codebook(:,symbol_indices_A1(j));
%             SCMA_codewords_A2(:,j)=present_codebook(:,symbol_indices_A2(j));
%             SCMA_codewords_A3(:,j)=present_codebook(:,symbol_indices_A3(j));
%             SCMA_codewords_A4(:,j)=present_codebook(:,symbol_indices_A4(j));
%             
%             SCMA_codewords_A5(:,j)=present_codebook(:,symbol_indices_A5(j));
%             SCMA_codewords_A6(:,j)=present_codebook(:,symbol_indices_A6(j));
%             SCMA_codewords_A7(:,j)=present_codebook(:,symbol_indices_A7(j));
%             SCMA_codewords_A8(:,j)=present_codebook(:,symbol_indices_A8(j));
%             
        end
         SCMA_codewords=SCMA_codewords_A1   ;
        
%            SCMA_codewords=[SCMA_codewords_A1 SCMA_codewords_A2 SCMA_codewords_A3 SCMA_codewords_A4...
%                            SCMA_codewords_A5 SCMA_codewords_A6 SCMA_codewords_A7 SCMA_codewords_A8];
            SCMA_codewords_MIMO=zeros(Nr*K,Nt*J);        
        for nr=1:Nr
            SCMA_codewords_MIMO((nr-1)*K+1:nr*K,:)=SCMA_codewords;
        end 
        symbol_indices=symbol_indices_A1 ;
%         symbol_indices=[symbol_indices_A1 symbol_indices_A2...
%                         symbol_indices_A3 symbol_indices_A4...
%                         symbol_indices_A5 symbol_indices_A6...
%                         symbol_indices_A7 symbol_indices_A8];
%   
        
        H_matrix=1/sqrt(2)*(randn(Nr*K,Nt*J)+1j*randn(Nr*K,Nt*J));
%         for nr=1:Nr
%         SCMA_codewords_MA((nr-1)*K+1:nr*K,:)=SCMA_codewords;
%         H_matrix((nr-1)*K+1:nr*K,:)=1/sqrt(2)*(randn(K,J)+1j*randn(K,J));
%         end
        %% Transmission through Rayleigh fading channel %%
        AWGN=sigma(e)*(randn(K*Nr,1)+1j*randn(K*Nr,1));    % complex Gaussian noise
        
       % h_matrix=1/sqrt(2)*(randn(K,J)+1j*randn(K,J));  % complex Rayleigh fading coefficient vector for all users
        %h_matrix=ones(K,J);   % Simple AWGN channel
        %h=1/sqrt(2)*(randn(K,1)+1j*randn(K,1));
       %h_matrix=repmat(h,1,J);% complex Rayleigh fading coefficient vector for Down-link        
        y=sum(H_matrix.*SCMA_codewords_MIMO,2)+AWGN;
        %y=reshape(y,[K,Nr]);
        %% received SCMA codeword UP_link         
       
        %[Rx_ant_prod,T1,T2,V,U,V_posterior,symbol_indices_hat]=SCMA_detection_MPA_prob_domain_MU_MIMO(power_users,y,N0(e),H_matrix,F,max_iter,M,C,Nr);
        %[U,V,arg1,T1,T2,V_posterior,symbol_indices_hat]=SCMA_detection_MPA_LLR_domain(power_users,y,N0(e),H_matrix,F1,max_iter,M,C,cf,Nr);
        [P,U,APP,symbol_indices_hat]=GAI_MIMO_SCMA_detection_real_form_mex(y,H_matrix,C,N0(e),J,K,Nt,Nr,M,T,F,F1);
        error_locations=find(symbol_indices~=symbol_indices_hat);
        
        %demapping
        %symbol_indices_hat=symbol_indices_hat-1;
%         bits_hat=de2bi(symbol_indices_hat,'left-msb');
%         bits_error=sum(sum(bits~=bits_hat));
        if ~isempty(error_locations)
            block_errors=block_errors+1;
            symbol_errors=symbol_errors+length(error_locations);
            %fprintf('\n  %d error collected',block_errors); 
        end     
        block=block+1;
    end
    SER(e)=symbol_errors/block/J/Nt;% Eack block contains J-symbols
    %BER(e)=SER(e)/m;
    fprintf('\n Simulation done for %d dB',Eb_N0_dB(e)); 
end
semilogy(Eb_N0_dB,SER,'b-*','LineWidth',2) ;
xlabel('Eb/No');
ylabel('SER');
legend('Nr=8,nt=1,M=4 AMPA');





