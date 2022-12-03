Purpose
======
We hope to study the functional relationship between adaptive phenomena and cortical networks in the primary visual cortex. For this reason, we have adopted a full model of our predecessors.
In the first step, we tried to simplify the model to get a large number of results quickly.

# full model
\> main_mutimodel  

# IVMM: Reduced Model 1
\> load('full_model_res.mat')  
\> \[] = compute_IVrefMM(res_full_model,params);  
\> load('full_model_res.mat2')  
\> \[estate_matref,istate_matref] = compute_IVrefMM(res_full_model,params, estate_matref,istate_matref);  
> Run IVMM (This step can be directly carried out when there is a transfer matrix “eirefstate_mat.mat”)  
\> IVrefMM_list = generate_IVrefMM(estate_matref,istate_matref);  
\> plot_fft(IVrefMM_list(1,:),IVrefMM_list(2,:));  
