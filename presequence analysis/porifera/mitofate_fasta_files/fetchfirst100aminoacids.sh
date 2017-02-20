for file in *.cleaned;
do
sed 's/*//g' $file> "$file".char
perl seqEnds.pl -l 50,1 "$file".char | sed 's/--[A-Z]//g' > "$file".first50.fasta
done
wait
rm *.char
echo "done!"
