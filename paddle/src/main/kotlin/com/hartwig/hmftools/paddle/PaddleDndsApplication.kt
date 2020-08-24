package com.hartwig.hmftools.paddle

import com.hartwig.hmftools.paddle.cohort.CohortLoad
import com.hartwig.hmftools.paddle.dnds.DndsCvGene
import com.hartwig.hmftools.paddle.dnds.DndsMutation
import com.hartwig.hmftools.paddle.likelihood.LikelihoodGene
import com.hartwig.hmftools.paddle.mutation.MutationsGene
import org.apache.logging.log4j.LogManager

fun main(args: Array<String>) {

    try {
        PaddleDndsApplication().use { x -> x.run() }
    } catch (e: Exception) {
        println(e)
    }
}

class PaddleDndsApplication : AutoCloseable, Runnable {

    companion object {
        val logger = LogManager.getLogger(this::class.java)
    }

    private val startTime = System.currentTimeMillis()


    override fun run() {
        val cohortFile = "/Users/jon/hmf/repos/hmftools/paddle/src/main/resources/HmfTMB.tsv"
        val dndsCVFile = "/Users/jon/hmf/repos/hmftools/paddle/src/main/resources/HmfRefCDSCv.tsv"
        val mutationsFile = "/Users/jon/hmf/repos/hmftools/paddle/src/main/resources/DndsMutations.tsv"

        logger.info("Loading dNdScv values: $dndsCVFile")
        val dndsCv = DndsCvGene.fromFile(dndsCVFile).associateBy { x -> x.gene }

        logger.info("Loading cohort: $cohortFile")
        val cohortLoad = CohortLoad.fromFile(cohortFile)

        logger.info("Loading mutations: $mutationsFile")
        val dndsMutations = DndsMutation.fromFile(mutationsFile)

        logger.info("Calculating gene mutation summary")
        val oncoGeneMutations = MutationsGene.oncoGeneMutations(dndsMutations).associateBy { x -> x.gene }
        val tsgGeneMutations = MutationsGene.tsgGeneMutations(dndsMutations).associateBy { x -> x.gene }

        val oncoLikelihood = LikelihoodGene(cohortLoad, dndsCv, oncoGeneMutations)
        val tsgLikelihood = LikelihoodGene(cohortLoad, dndsCv, tsgGeneMutations)

        LikelihoodGene.writeFile(false, "/Users/jon/hmf/repos/hmftools/hmf-common/src/main/resources/dnds/DndsDriverLikelihoodOnco.tsv", oncoLikelihood.values)
        LikelihoodGene.writeFile(false, "/Users/jon/hmf/repos/hmftools/hmf-common/src/main/resources/dnds/DndsDriverLikelihoodTsg.tsv", tsgLikelihood.values)
    }



    override fun close() {
        logger.info("Finished in ${(System.currentTimeMillis() - startTime) / 1000} seconds")
    }
}