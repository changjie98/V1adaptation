Purpose
======
We hope to study the functional relationship between adaptive phenomena and cortical networks in the primary visual cortex. For this reason, we have adopted a full model of our predecessors.
In the first step, we tried to simplify the model to get a large number of results quickly.

## full model
\> main_mutimodel  

## IVMM: Reduced Model 1
### method
Replace variance information with mean value information of current and voltage  

### code useage
> Calculate the probability transfer matrix  
\> load('full_model_res.mat')  
\> \[estate_mat,istate_mat] = compute_IVrefMM(res_full_model,params);  
\> load('full_model_res.mat2')  
\> \[estate_mat,istate_mat] = compute_IVrefMM(res_full_model,params,estate_mat,istate_mat);  
> Run IVMM (This step can be directly carried out when there is a transfer matrix “eirefstate_mat.mat”)  
\> IVrefMM_list = generate_IVrefMM(estate_mat,istate_mat);  
\> plot_fft(IVrefMM_list(1,:),IVrefMM_list(2,:));  

## DIF-ODE: Reduced Model 2
### method
Mean field and ODE. Treat neurons with similar voltages as a distribution

### code useage
> Run the model
\> main_mutimodel
> Calculate Jacobian matrix
\> \[J_mat,D_vac] = get_jacobian_matrix(res_ODEfull_model,params, time_index);
