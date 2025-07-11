for r_id in 1 2 3 4 5 6 7 8 9 10 
do	
	for f_id in 1 2 3 4 5 6 7 8 9 10
	do
		for on_what in testing
			do
			for which_model in RF
				do
				for dataset in EO-CRC LO-CRC Alldata
				do
					echo -e "Rscript run_SHAP_test.R ${r_id} ${f_id} ${on_what} ${which_model} ${dataset}"
				done
			done
		done
	done
done > SHAP_job_list.sh
