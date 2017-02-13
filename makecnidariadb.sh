for file in *.transdecoder.pep;
do
grep ">" $file | sed 's/>//g' | awk -F "ORF" '{print $1}' > transdecoder.list 
done

for file in *.transdecoder_renamed.pep;
do
grep ">" $file | sed 's/>//g' > transdecoder.renamed.list 
done

paste -d "," transdecoder.list transdecoder.renamed.list > transcript_protein_renamed_ID_mapping.list

wait

while read p; do

for file in *.targetP.clean;
do
grep -w $p $file | awk '{print $6,",",$7}' >> targetP.map
done

done < transdecoder.renamed.list

wait

for file in *.human;
do
awk '{print $1}' $file | sort > "$file".humanlist.sorted
sort transdecoder.renamed.list > transdecoder.renamed.list.sorted
comm -23 transdecoder.renamed.list.sorted "$file".humanlist.sorted > notpresentinhuman.list
awk '{print $1,"\t","-"}' notpresentinhuman.list > notpresentinhuman.list.updated
cat notpresentinhuman.list.updated $file > FINAL.HUMAN  
done


for file in *.mouse;
do
awk '{print $1}' $file | sort > "$file".mouselist.sorted
sort transdecoder.renamed.list > transdecoder.renamed.list.sorted
comm -23 transdecoder.renamed.list.sorted "$file".mouselist.sorted > notpresentinmouse.list
awk '{print $1,"\t","-"}' notpresentinmouse.list > notpresentinmouse.list.updated
cat notpresentinmouse.list.updated $file > FINAL.MOUSE  
done

for file in *.yeast;
do
awk '{print $1}' $file | sort > "$file".yeastlist.sorted
sort transdecoder.renamed.list > transdecoder.renamed.list.sorted
comm -23 transdecoder.renamed.list.sorted "$file".yeastlist.sorted > notpresentinyeast.list
awk '{print $1,"\t","-"}' notpresentinyeast.list > notpresentinyeast.list.updated
cat notpresentinyeast.list.updated $file > FINAL.YEAST  
done


wait

while read p; do

for file in *FINAL.HUMAN;
do
grep -w $p $file | awk '{print $1,",",$2}' >> "$file".RBBH.map
done

done < transdecoder.renamed.list

wait

while read p; do

for file in *FINAL.MOUSE;
do
grep -w $p $file | awk '{print $1,",",$2}' >> "$file".RBBH.map
done

done < transdecoder.renamed.list

wait

while read p; do

for file in *FINAL.YEAST;
do
grep -w $p $file | awk '{print $1,",",$2}' >> "$file".RBBH.map
done

done < transdecoder.renamed.list

wait


awk -F "," '{print $2}' FINAL.HUMAN.RBBH.map > human.RBBH.results.final.map
awk -F "," '{print $2}' FINAL.MOUSE.RBBH.map > mouse.RBBH.results.final.map
awk -F "," '{print $2}' FINAL.YEAST.RBBH.map > yeast.RBBH.results.final.map

paste -d "," transcript_protein_renamed_ID_mapping.list targetP.map human.RBBH.results.final.map mouse.RBBH.results.final.map yeast.RBBH.results.final.map > IDmap_TargetPmap_RBBHmap.list

for file in *uni_reciprocal_targetP.LST.prot;
do
grep ">" $file | sed 's/>//g' | sort > mitoprot.list.sorted
done

wait

comm -23 transdecoder.renamed.list.sorted mitoprot.list.sorted | sort > nonmitprot.list.sorted

wait

awk '{print $1,",","Y"}' mitoprot.list.sorted > mitoprot.updated
awk '{print $1,",","N"}' nonmitprot.list.sorted > nonmitoprot.updated
cat mitoprot.updated nonmitoprot.updated > mitoprot.updated.map


wait

while read p; do

for file in *mitoprot.updated.map;
do
grep -w $p $file | awk -F "," '{print $1,",",$2}' >> "$file".mitoprot.map.final
done

done < transdecoder.renamed.list

awk -F "," '{print $2}' mitoprot.updated.map.mitoprot.map.final > mitoprot.map

paste -d "," IDmap_TargetPmap_RBBHmap.list mitoprot.map > IDmap_TargetPmap_RBBHmap_Mitoprotmap.list 


grep ">" Human.MitoCarta2.0.fasta | sed 's/>//g' | awk '{print $3}' > humanmitocarta.genelist
grep ">" human_renamed.pep | sed 's/>//g' > humanmitocarta.renamed.list
paste -d "," humanmitocarta.renamed.list humanmitocarta.genelist > human_mitocarta_gene.map

for file in *.human;
do
awk '{print $2}' $file > "$file".humanlist.forgene


while read p;do
grep -w $p human_mitocarta_gene.map | awk -F "," '{print $2}' >> genemap
done < "$file".humanlist.forgene

done

wait

for file in *.human;
do
paste -d "\t" $file genemap > human.genemap.updated
done


for file in *.human;
do
awk '{print $1}' $file | sort > "$file".humanlist.sorted.forgene
sort transdecoder.renamed.list > transdecoder.renamed.list.sorted.forgene
comm -23 transdecoder.renamed.list.sorted.forgene "$file".humanlist.sorted.forgene > notpresentinhuman.list.forgene
awk '{print $1,"\t","-","\t","-"}' notpresentinhuman.list.forgene > notpresentinhuman.list.updated.forgene
cat notpresentinhuman.list.updated.forgene human.genemap.updated > FINAL.HUMAN.forgene  
done


while read p; do

for file in *FINAL.HUMAN.forgene;
do
grep -w $p $file | awk '{print $1,",",$2,","$3}' >> "$file".RBBH.map.forgene
done

done < transdecoder.renamed.list


awk -F "," '{print $3}' FINAL.HUMAN.forgene.RBBH.map.forgene > genemitocarta.map

paste -d "," IDmap_TargetPmap_RBBHmap_Mitoprotmap.list genemitocarta.map > IDmap_TargetPmap_RBBHmap_Mitoprotmap_Genemap.list 

for f in *transdecoder_renamed.pep;
do
perl fasta_get_seq_length.pl < $f > "$f".seq-length
paste -d "," IDmap_TargetPmap_RBBHmap_Mitoprotmap_Genemap.list "$f".seq-length > IDmap_TargetPmap_RBBHmap_Mitoprotmap_Genemap_Seqlength.list

done

wait
grep ">" MouseMitoCarta.proteins.fasta | sed 's/>//g' | awk '{print $3}' > Mousemitocarta.genelist
grep ">" mouse_renamed.pep | sed 's/>//g' > Mousemitocarta.renamed.list
paste -d "," Mousemitocarta.renamed.list Mousemitocarta.genelist > mouse_mitocarta_gene.map


for file in *.mouse;
do
awk '{print $2}' $file > "$file".mouselist.forgene


while read p;do
grep -w $p mouse_mitocarta_gene.map | awk -F "," '{print $2}' >> mouse.genemap
done < "$file".mouselist.forgene

done

wait

for file in *.mouse;
do
paste -d "\t" $file mouse.genemap > mouse.genemap.updated
done

wait


for file in *.mouse;
do
awk '{print $1}' $file | sort > "$file".mouselist.sorted.forgene
sort transdecoder.renamed.list > transdecoder.renamed.list.sorted.forgene
comm -23 transdecoder.renamed.list.sorted.forgene "$file".mouselist.sorted.forgene > notpresentinmouse.list.forgene
awk '{print $1,"\t","-","\t","-"}' notpresentinmouse.list.forgene > notpresentinmouse.list.updated.forgene
cat notpresentinmouse.list.updated.forgene mouse.genemap.updated > FINAL.mouse.forgene  
done

wait

while read p; do

for file in *FINAL.mouse.forgene;
do
grep -w $p $file | awk '{print $1,",",$2,","$3}' >> "$file".RBBH.map.forgene
done

done < transdecoder.renamed.list



awk -F "," '{print $3}' FINAL.mouse.forgene.RBBH.map.forgene > genemitocarta.mouse.map

paste -d "," IDmap_TargetPmap_RBBHmap_Mitoprotmap_Genemap_Seqlength.list genemitocarta.mouse.map > IDmap_TargetPmap_RBBHmap_Mitoprotmap_Genemap_Seqlength_MouseGenemap.list


grep ">" yeast_mitochondrial.fasta | awk '{print $2}' > yeastmitochondria.genelist
grep ">" yeast_mitochondrial_renamed.fasta | sed 's/>//g' > yeast_mitochondrial.renamed.list
paste -d "," yeast_mitochondrial.renamed.list yeastmitochondria.genelist > yeast_gene.map


for file in *.yeast;
do
awk '{print $2}' $file > "$file".yeastlist.forgene


while read p;do
grep -w $p yeast_gene.map | awk -F "," '{print $2}' >> yeast.genemap
done < "$file".yeastlist.forgene

done

wait

for file in *.yeast;
do
paste -d "\t" $file yeast.genemap > yeast.genemap.updated
done

wait


for file in *.yeast;
do
awk '{print $1}' $file | sort > "$file".yeastlist.sorted.forgene
sort transdecoder.renamed.list > transdecoder.renamed.list.sorted.forgene
comm -23 transdecoder.renamed.list.sorted.forgene "$file".yeastlist.sorted.forgene > notpresentinyeast.list.forgene
awk '{print $1,"\t","-","\t","-"}' notpresentinyeast.list.forgene > notpresentinyeast.list.updated.forgene
cat notpresentinyeast.list.updated.forgene yeast.genemap.updated > FINAL.yeast.forgene  
done

wait

while read p; do

for file in *FINAL.yeast.forgene;
do
grep -w $p $file | awk '{print $1,",",$2,","$3}' >> "$file".RBBH.map.forgene
done

done < transdecoder.renamed.list


awk -F "," '{print $3}' FINAL.yeast.forgene.RBBH.map.forgene > genemitocarta.yeast.map

paste -d "," IDmap_TargetPmap_RBBHmap_Mitoprotmap_Genemap_Seqlength_MouseGenemap.list genemitocarta.yeast.map > IDmap_TargetPmap_RBBHmap_Mitoprotmap_Genemap_Seqlength_MouseGenemap_YeastGenemap.list
