Purpose
======
We hope to study the functional relationship between adaptive phenomena and cortical networks in the primary visual cortex. For this reason, we have adopted a full model of our predecessors.
In the first step, we tried to simplify the model to get a large number of results quickly.

# full model
/> main_mutimodel

# IVMM: Reduced Model 1
&ltload('full_model_res.mat')
&lt[] = compute_IVrefMM(res_full_model,params);
&ltload('full_model_res.mat2')
&lt[estate_matref,istate_matref] = compute_IVrefMM(res_full_model,params, estate_matref,istate_matref);
Run IVMM (This step can be directly carried out when there is a transfer matrix “eirefstate_mat.mat”)
&ltIVrefMM_list = generate_IVrefMM(estate_matref,istate_matref);
&ltplot_fft(IVrefMM_list(1,:),IVrefMM_list(2,:));
