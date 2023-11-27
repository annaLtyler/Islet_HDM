## Run multiple gene sets through 4.Gene_Set_Validation.Rmd
## and rename the results 



##=====================================================================================##
## Leptin
##=====================================================================================##

R -e "rmarkdown::render(here::here('Documents', '4.Gene_Set_Validation.Rmd'))" --args "Lep-Lepr_quaternary_protein_cluster.txt"
mv Documents/4.Gene_Set_Validation.html Documents/4.Gene_Set_Validation_Leptin.html


##=====================================================================================##
## Apelin
##=====================================================================================##

R -e "rmarkdown::render(here::here('Documents', '4.Gene_Set_Validation.Rmd'))" --args "Apelin_signaling_pathway_cluster.txt"
mv Documents/4.Gene_Set_Validation.html Documents/4.Gene_Set_Validation_Apelin.html


##=====================================================================================##
## endomembrane system
##=====================================================================================##

R -e "rmarkdown::render(here::here('Documents', '4.Gene_Set_Validation.Rmd'))" --args "endomembrane_system_cluster.txt"
mv Documents/4.Gene_Set_Validation.html Documents/4.Gene_Set_Validation_Endomembrane_System.html

