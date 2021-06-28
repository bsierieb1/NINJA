# NINJA (Nematode INdicator Joint Analysis) version 1

This is a complete Shiny app which may be run locally (check out the tutorials at https://shiny.rstudio.com/tutorial/ if you are not familiar with Shiny). Make sure that ui.R, server.R, database.xls, and www are in the same directory.

If you only need to access certain fragments of the code (e.g. a function to calculate the Maturity Index), look in the file server.R.

In contrast to the NINJA code, the database containing all the information about taxa is a reformatted version of the Nemaplex database (http://nemaplex.ucdavis.edu/). It is intellectual property of Howard Ferris and we have no permission to share it. Therefore, you have to create your own database file. You may email us for tips. The database.xls file in this repository is a demo version which we added as an example of correct formatting. The columns "taxon", "cp", "feeding", and "mass" are self-explanatory. As for the other columns, "plant" is a subcategory of the feeding type 1 according to Yeates et al. (i.e. for feeding type 1b one should type "1" in the column "feeding" and "b" in the column "plant"). For non-plant feeders, this column should be left empty. "level" is taxonomic level, where "f" means family, "g" - genus, and "s" - species. The columns "fung_family" and "diplogasterid" are only relevant if you need to calculate the compost maturity index. If you do not use this metric, you may as well leave these fields blank. Please make sure that you keep the taxa whose names start with "zzzfake..." in the database file!

Please cite: https://doi.org/10.1016/j.ejsobi.2014.02.004

Contact: sarasm@inia.es
