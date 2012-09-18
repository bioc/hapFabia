#
#
# Author: SEPP HOCHREITER
###############################################################################

setClass("HaploCluster",
         representation = representation(
         ID = "numeric",
         bicluster_id = "numeric",
         chromosome = "character",
         haploClusterPos = "numeric",
         haploClusterLength = "numeric",
         numberIndividuals = "numeric",
         numbertagSNVs = "numeric",
         individuals = "vector",
         tagSNVs = "vector",
         populationIndividuals = "vector",
         idIndividuals = "vector",
         labelIndividuals = "vector",
         platformIndividuals = "vector",
         coreClusterIndividuals = "vector",
         tagSNVPositions = "vector",
         tagSNVAlleles = "vector",
         tagSNVNames = "vector",
         tagSNVFreq = "vector",
         tagSNVGroupFreq = "vector",
         tagSNVChange = "vector",
         tagSNVsPerIndividual = "vector",
         individualPerTagSNV = "vector",
         tagSNVAnno = "vector"
         )
)


setValidity("HaploCluster",
    function(object)
    {
        if (!is.numeric(slot(object, "ID")) || (slot(object, "ID")<0))
        {
            return("slot >ID< must be an integer number larger 0!")
        }
        else if (!is.numeric(slot(object, "bicluster_id")) || (slot(object, "bicluster_id")<0))
        {
            return("slot >bicluster_id< must be an integer number larger 0!")
        }
        else if (!is.character(slot(object, "chromosome")))
        {
            return("slot >chromosome< must be a string!")
        }
        else if (!is.numeric(slot(object, "haploClusterPos")) || (slot(object, "haploClusterPos")<0))
        {
            return("slot >haploClusterPos< must be an integer number larger 0!")
        }
        else if (!is.numeric(slot(object, "haploClusterLength")) || (slot(object, "haploClusterLength")<0))
        {
            return("slot >haploClusterLength< must be an integer number larger 0!")
        }
        else if (!is.numeric(slot(object, "numberIndividuals")) || (slot(object, "numberIndividuals")<0))
        {
            return("slot >numberIndividuals< must be an integer number larger 0!")
        }
        else if (!is.numeric(slot(object, "numbertagSNVs")) || (slot(object, "numbertagSNVs")<0))
        {
            return("slot >numbertagSNVs< must be an integer number larger 0!")
        }
        else if (!is.vector(slot(object, "individuals")))
        {
            return("slot >individuals< must be a vector!")
        }
        else if (!is.vector(slot(object, "tagSNVs")))
        {
            return("slot >tagSNVs< must be a vector!")
        }
        else if (!is.vector(slot(object, "populationIndividuals")))
        {
            return("slot >populationIndividuals< must be a character vector!")
        }
        else if (!is.vector(slot(object, "idIndividuals")))
        {
            return("slot >idIndividuals< must be a vector!")
        }
        else if (!is.vector(slot(object, "labelIndividuals")))
        {
            return("slot >labelIndividuals< must be a character vector!")
        }
        else if (!is.vector(slot(object, "platformIndividuals")))
        {
            return("slot >platformIndividuals< must be a character vector!")
        }
        else if (!is.vector(slot(object, "coreClusterIndividuals")))
        {
            return("slot >coreClusterIndividuals< must be a vector!")
        }
        else if (!is.vector(slot(object, "tagSNVPositions")))
        {
            return("slot >tagSNVPositions< must be a vector!")
        }
        else if (!is.vector(slot(object, "tagSNVAlleles")))
        {
            return("slot >tagSNVAlleles< must be a character vector!")
        }
        else if (!is.vector(slot(object, "tagSNVNames")))
        {
            return("slot >tagSNVNames< must be a character vector!")
        }
        else if (!is.vector(slot(object, "tagSNVFreq")))
        {
            return("slot >tagSNVFreq< must be a vector!")
        }
        else if (!is.vector(slot(object, "tagSNVGroupFreq")))
        {
            return("slot >tagSNVGroupFreq< must be a vector!")
        }
        else if (!is.vector(slot(object, "tagSNVChange")))
        {
            return("slot >tagSNVChange< must be a vector!")
        }
        else if (!is.vector(slot(object, "tagSNVsPerIndividual")))
        {
            return("slot >tagSNVsPerIndividual< must be a vector!")
        }
        else if (!is.vector(slot(object, "individualPerTagSNV")))
        {
            return("slot >individualPerTagSNV< must be a vector!")
        }
        else if (!is.vector(slot(object, "tagSNVAnno")))
        {
            return("slot >tagSNVAnno< must be a character vector!")
        }
    }
 )


setClass("HaploClusterList",
         representation = representation(
         haploClusters="list",
         lengthList="numeric",
         statistics="list"
         )
)


setValidity("HaploClusterList",
    function(object)
    {
        if (!is.list(slot(object, "haploClusters")))
        {
            return("slot >haploClusters< must be a list!")
        } else if (!is.numeric(slot(object, "lengthList")) || (slot(object, "lengthList")<0))
        {
            return("slot >lengthList< must be an integer number larger equal 0!")
        }  else if (!is.list(slot(object, "statistics")))
        {
            return("slot >statistics< must be a list!")
        }


    }

 )

