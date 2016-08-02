for date in data/*
do 
	for seq in $date/*
	do
		if [ ! -d $seq ]; then
			continue
		fi	
			
		for evtfile in $seq/event_cl/*N_cl_sunpos.evt
		do
			./nustar_make_sunpy_map.py -i $evtfile
		done
		

	done	
done
