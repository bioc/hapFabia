# Init file for package fabia
.onLoad <- function(lib, pkg)
{

    library.dynam("hapFabia", pkg, lib)

    if ((.Platform$OS.type == "windows") && (.Platform$GUI ==
        "Rgui") && interactive()) {
        vigFile = system.file("Meta", "vignette.rds", package = "hapFabia")
        if (!file.exists(vigFile)) {
            warning(sprintf("hapFabia vignette is missing, nothing is added to the menu bar"))
        }
        else {
            vigMtrx = readRDS(vigFile)
            vigs = file.path(chartr("\\", "/", find.package("hapFabia")), "doc", vigMtrx[,
                "PDF"])
            names(vigs) = vigMtrx[, "Title"]
            if (!"Vignettes" %in% winMenuNames())
                winMenuAdd("Vignettes")
            pkgMenu = paste("Vignettes", "hapFabia", sep = "/")
            winMenuAdd(pkgMenu)
            for (i in seq(along = vigs)) winMenuAddItem(pkgMenu,
                names(vigs)[i], paste("shell.exec(\"", vigs[i],
                  "\")", sep = ""))
        }
    }

  packageStartupMessage(
      "+--------------------------+          #    #    ##    #####     \n",
      "|#.....#...#.......#.#....#|          #    #   #  #   #    #    \n",
      "|#.....#...#.......#.#....#|          ######  #    #  #    #    \n",
      "|#.....#...#...............|          #    #  ######  #####     \n",
      "|#.....#...#.......#.#....#|          #    #  #    #  #         \n",
      "|#.....#...#...............|          #    #  #    #  #         \n",
      "|#.....#...#.......#.#....#|  #######                           \n",
      "|..................#.#....#|  #         ##    #####   #    ##   \n",
      "|#.....#...#.......#.#....#|  #        #  #   #    #  #   #  #  \n",
      "|..................#.#....#|  #####   #    #  #####   #  #    # \n",
      "|#.....#...#.......#.#....#|  #       ######  #    #  #  ###### \n",
      "|#.....#...#.......#.#....#|  #       #    #  #    #  #  #    # \n",
      "+--------------------------+  #       #    #  #####   #  #    # \n")
    version <- packageDescription("hapFabia",fields="Version")
    packageStartupMessage( "Citation: S. Hochreiter,","\n",
      "HapFABIA: Identification of very short segments of identity by descent characterized by rare variants in large sequencing data,","\n",
      "Nucleic Acids Research, 2013, doi: 10.1093/nar/gkt1013.","\n",
      "BibTex: enter 'toBibtex(citation(\"hapFabia\"))'","\n\n",
      "Homepage: http://www.bioinf.jku.at/software/hapFabia/index.html","\n\n",
      "hapFabia Package Version ", version, "\n")
}

.onUnload <- function(libpath)
{
    library.dynam.unload("hapFabia", libpath)
}


