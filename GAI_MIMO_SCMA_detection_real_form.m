function [P,U,APP,symbol_indices_hat]=GAI_MIMO_SCMA_detection_real_form(y,H,C,sigma_v,J,K,Nt,Nr,M,T,F,F1)
P_in=repmat(F1',M,1); 
P_in=reshape(P_in,J*Nt,M*K*Nr);
P_in=1/M*P_in';
P=P_in;
U=zeros(K*M*Nr,J*Nt);
%V=zeros(K*M,J);
%mean_mat=zeros(Nr*K,Nt*J);
%var_mat=zeros(Nr*K,Nt*J);
APP=zeros(M,J*Nt);
for t=1:T
    APP_prev=APP;   
 %Compu of APP at FN
 for nr=1:Nr
     H_nr=H((nr-1)*K+1:nr*K,:);%present nr channel matrix
     y_nr=y((nr-1)*K+1:nr*K);%present nr obser vector
   for k=1:K
       y_nr_k=y_nr(k);%k-th observation 
       y_nr_k_RV=[real(y_nr_k);imag(y_nr_k)];%real form
       NE_VN=find(F(k,:)); %neigh VNs
       S_mean_k=zeros(2,1);
       K_matrix=zeros(2,2);
      for nt=1:Nt 
            H_nr_nt=H_nr(:,(nt-1)*J+1:nt*J);%present (nr,nt) channel matrix
            pr_nr_nt_P=P((nr-1)*K*M+1:nr*K*M,(nt-1)*J+1:nt*J);%present nt and nr codeword Probablities            
       for j=1:length(NE_VN)
            pr_j=NE_VN(j);
            pr_CB=C((pr_j-1)*K+1:pr_j*K,:);%pr_CB
            pr_cw_RE=pr_CB(k,:);%pr_RE            
            pr_k_P=pr_nr_nt_P((k-1)*M+1:k*M,pr_j);%present RE priori
            mean_k=zeros(2,1);
            K_pr_j=zeros(2,2);
            for m=1:M
            v_j_k=H_nr_nt(k,pr_j)*pr_cw_RE(m);%complex scalar
            Real_v_j_k=[real(v_j_k);imag(v_j_k)];%real vector form
            mat=Real_v_j_k*Real_v_j_k';%real matrix for computing covariance
            mean_k=mean_k+pr_k_P(m)*Real_v_j_k;
            K_pr_j= K_pr_j+pr_k_P(m)*mat;
            end
            K_pr_j=K_pr_j-mean_k*mean_k';%pr_UE Covariance
            %mean_mat((k-1)*2+1:k*2,pr_j)=mean_k;
            S_mean_k=S_mean_k+mean_k;
            K_matrix=K_matrix+K_pr_j;                     
       end
      end 
   for nt=1:Nt
       H_nr_nt=H_nr(:,(nt-1)*J+1:nt*J);%present (nr,nt) channel matrix
       pr_nr_nt_P=P((nr-1)*K*M+1:nr*K*M,(nt-1)*J+1:nt*J);
       for j=1:length(NE_VN)
            pr_j=NE_VN(j);%present UE on k-th RE
            pr_CB=C((pr_j-1)*K+1:pr_j*K,:);%pr_CB
            pr_cw_RE=pr_CB(k,:);%pr_RE
            pr_k_P=pr_nr_nt_P((k-1)*M+1:k*M,pr_j);
            pr_mean_k=zeros(2,1);
            pr_K=zeros(2,2);
            for m=1:M
            v_j_k=H_nr_nt(k,pr_j)*pr_cw_RE(m);
            Real_v_j_k=[real(v_j_k);imag(v_j_k)];%real vector form
            mat=Real_v_j_k*Real_v_j_k';%real matrix for cov_cal
            pr_mean_k=pr_mean_k+pr_k_P(m)*Real_v_j_k;%pr_mean
            pr_K=pr_K+pr_k_P(m)*mat;%pr_covariance
            end
            pr_K=pr_K-pr_mean_k*pr_mean_k';
            S_mu=S_mean_k-pr_mean_k;%mean vector of (I+N)
            S_K=K_matrix-pr_K+sigma_v/2*eye(2,2);%covariance matrix of (I+N)
            %S_K=S_K+0.3*eye(2);
            U_message=zeros(M,1);
            for m=1:M
             %pr_CW=pr_CB(:,m);
            v_j_k=H_nr_nt(k,pr_j)*pr_cw_RE(m);
            Real_v_j_k=[real(v_j_k);imag(v_j_k)];
            U_message(m)=(1/(2*pi*sqrt(det(S_K))))*exp((-1/2)*(y_nr_k_RV-Real_v_j_k-S_mu)'*(eye(2)/S_K)*(y_nr_k_RV-Real_v_j_k-S_mu));
            % U_message(m)=exp((-1/2)*(y_nr_k_RV-Real_v_j_k-S_mu)'*(eye(2)/S_K)*(y_nr_k_RV-Real_v_j_k-S_mu));
             
            end
             U_message=U_message/sum(U_message);
             U((nr-1)*M*K+(k-1)*M+1:(nr-1)*K*M+k*M,(nt-1)*J+NE_VN(j))=U_message;
            
        end  
   end     
   end
 end
 
 %VN update
 for nt=1:Nt
 for j=1:J
     NE_RN=find(F(:,j)); %Neigh RNs  
    for nr=1:Nr% present nr

     for k=1:length(NE_RN)% present k
         pr_k=NE_RN(k);
         NE_RN_min_k=setdiff(NE_RN,pr_k);
         P_message=zeros(M,1);         
         for m=1:M
             p1=1;
            for nr1=1:Nr
                 U_nr_nt=U((nr1-1)*K*M+1:nr1*K*M,(nt-1)*J+1:nt*J);
             for i1=1:length(NE_RN_min_k)
                 pr_k1=NE_RN_min_k(i1);
                 U_nr_pr_k=U_nr_nt((pr_k1-1)*M+1:pr_k1*M,j);
                 p1=p1*U_nr_pr_k(m);                
             end
            end
             P_message(m)=p1;  
            % P((nr-1)*K*M+(pr_k-1)*M+m,(nt-1)*J+j)=p1;
             
         end  
         P_message=P_message/sum(P_message);
         P((nr-1)*K*M+(pr_k-1)*M+1:(nr-1)*M*K+pr_k*M,(nt-1)*J+j)=P_message;
         
     end
    end
 end
 end
 
 for nt=1:Nt
  for j=1:J
     NE_RN=find(F(:,j)); %Neigh RNs  
         APP_message=zeros(M,1);
         for m=1:M
             p1=1;
            for nr1=1:Nr
                U_nr_nt=U((nr1-1)*K*M+1:nr1*K*M,(nt-1)*J+1:nt*J);
             for i1=1:length(NE_RN)
                 pr_k1=NE_RN(i1);
                 U_nr_pr_k=U_nr_nt((pr_k1-1)*M+1:pr_k1*M,j);
                 p1=p1*U_nr_pr_k(m);                
             end
            end
             APP_message(m)=p1;  
             %APP(m,(nt-1)*J+j)=p1;
             
         end
         APP_message=APP_message/sum(APP_message);
         APP(:,(nt-1)*J+j)=APP_message;
         
         
  
  end
 end
 %%
 % Check for convergence
    check_value=sum(sum(abs(APP_prev-APP)));
    if check_value<0.0001         % stopping rule        
        break;
    end
end
%% Estimates of the codeword
symbol_indices_hat=zeros(J,Nt);
max_values=max(APP);
for nt=1:Nt
for j=1:J
    for m=1:M
        if APP(m,(nt-1)*J+j)==max_values((nt-1)*J+j)
            symbol_indices_hat(j,nt)=m;
            break;
        end
    end
end  
end