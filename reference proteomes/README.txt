This folder contains the original reference mitochondrial proteomes and the files with renamed headers.
Headers were renamed using 
awk '/^>/{print ">shys" ++i; next}{print}' <   file    > renamed.file
