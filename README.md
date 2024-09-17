Laura Cooke Scannells Master's Thesis: "Life Between a Rock and an Interface: a Tale of Two Caves"

This github page contains several examples of the R code I used for my master's thesis, which compares the diversity and ecology of two novel microbiological environments: Grotta Del Cane in Italy and Sulfur Cave in Romania. Both caves are examples of a little-studied ecosystem known as a gas-gas redox interface ecosystem. Gas-gas redox interface ecosystems are environments where two distinct layers of gasses form an interface between each other, which contains gasses from both layers and their redox interactions, for example, hydrogen sulfide and oxygen form a redox couple, which can be taken advantage of by microorganisms for energy.

Both caves have two layers of gas: a lower layer comprised of extremely high concentrations of carbon dioxide, and an upper layer of atmospheric air. The caves are extremely dry, limiting their ability to utilize water for movement and metabolism. Both caves are extremely acidic (pH < 1) and completely dark, forcing the organisms to live out their lives through chemosynthesis in extreme conditions.

To compare the microbial communities in both caves, my thesis team and I collected samples from both caves from the top of the cave walls to the bottom of the cave walls to compare the levels of stratification between both caves and characterize the microbial communities. Then I extracted the 16S NGS Illumina DNA from the samples and utilized a mix of Qiime2, R, and Excel to analyze them.

The "Qiime2 rundown" document is an example of the Qiime2 techniques I used for my 16S bioinformatics pipeline, which took raw sequencing data and converted it to .qza files that can be imported into R through the Qiime2r package (https://github.com/jbisanz/qiime2R) and then interpreted through familiar R packages such as Bioconductor, vegan, etc. Furthermore, the Qiime2 rundown also shows the use of Qiime2's own plugins, which can also be used to generate figures within Qimme2 itself such as the "alpha-group-significance" command.

The "Joostpoopunifrac" and "grottaunifrac" files are examples of the use of Qiime2 data in R to produce unifrac 
