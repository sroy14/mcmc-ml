aws_hmc1:
	o2graph -set fig_dict "dpi=250" \
		-read ../../o2scl/examples/data/ex_mcmc_nn1_0_out \
		-select-rows "class>0.5" -scatter x y \
		-ttext 0.3 0.9 "HMC samples" "ha=center" \
		-read ../../o2scl/examples/data/ex_mcmc_nn1_0_out \
		-select-rows "class<0.5" -scatter x y \
		-xtitle x -ytitle y -save aws_hmc1.pdf -show

aws_hmc2:
	o2graph -set fig_dict "dpi=250" \
		-read ../../o2scl/examples/data/ex_mcmc_nn1_0_out \
		-select-rows "class<0.5" -scatter x y \
		-ttext 0.3 0.9 "HMC class failures" "ha=center" \
		-xtitle x -ytitle y -save aws_hmc2.pdf -show

aws_prob1:
	o2graph -set fig_dict "dpi=250" \
		-subplots 2 2 \
		-subadj "top=0.99,right=0.99,left=0.14,hspace=0.28,wspace=0.35" \
		-selax 0 \
		-read ../../o2scl/examples/data/ex_mcmc_nn2.o2 prop_GMM \
		-scatter lwp lwq -xtitle exact -ytitle GMM \
		-set font 12 \
		-ttext 0.12 0.95 "(false positives)" "va=top,rotation=90" \
		-ttext 0.95 0.12 "(false negatives)" "ha=right" \
		-set font 16 \
		-selax 1 \
		-read ../../o2scl/examples/data/ex_mcmc_nn2.o2 prop_KDE \
		-scatter lwp lwq -xtitle exact -ytitle KDE \
		-selax 2 \
		-read ../../o2scl/examples/data/ex_mcmc_nn2.o2 prop_nflow \
		-scatter lwp lwq -xtitle exact -ytitle "norm. flow" \
		-selax 3 -eval "self.axes.axis('off')" \
		-save aws_prob1.pdf -show

aws_prob2:
	o2graph -set fig_dict "dpi=250" \
		-subplots 2 2 \
		-subadj "top=0.99,right=0.99,left=0.14,hspace=0.28,wspace=0.35" \
		-set xlo "(-3)" -set xhi 0.2 \
		-set ylo "(-3)" -set yhi 0.2 \
		-selax 0 \
		-read ../../o2scl/examples/data/ex_mcmc_nn2.o2 prop_GMM \
		-scatter lwp lwq -xtitle exact -ytitle GMM \
		-selax 1 \
		-read ../../o2scl/examples/data/ex_mcmc_nn2.o2 prop_KDE \
		-scatter lwp lwq -xtitle exact -ytitle KDE \
		-selax 2 \
		-read ../../o2scl/examples/data/ex_mcmc_nn2.o2 prop_nflow \
		-scatter lwp lwq -xtitle exact -ytitle "norm. flow" \
		-selax 3 -eval "self.axes.axis('off')" \
		-save aws_prob2.pdf -show




