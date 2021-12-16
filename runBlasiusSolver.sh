while read Me; do
	echo -e "$Me\n" | python -W ignore compressibleBlasius.py
# done <Me_buildSmartGuessList.txt
done <Me_select.txt
# echo -e "10\n" | python -W ignore compressibleBlasius.py
# python post_processing_compressibleBlasius.py