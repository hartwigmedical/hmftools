package com.hartwig.hmftools.paddle

import com.hartwig.hmftools.paddle.cohort.TumorMutationalLoad
import com.hartwig.hmftools.paddle.cohort.TumorMutationalLoadSample
import com.hartwig.hmftools.paddle.dnds.DndsCvGene
import com.hartwig.hmftools.paddle.dnds.DndsMutation
import com.hartwig.hmftools.paddle.driver.LikelihoodGene
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
        val mutationsFile = "/Users/jon/hmf/repos/hmftools/paddle/src/main/resources/DndsMutations.AR.tsv"

        logger.info("Loading cohort: $cohortFile")
        val (cohortSize, cohortLoad) = loadCohort(cohortFile)

        logger.info("Loading mutations: $mutationsFile")
//        val dndsMutations = DndsMutation.fromFile(mutationsFile)
        val dndsMutations = DndsMutation.fromFile(mutationsFile).filter { x -> x.impact == Impact.MISSENSE || x.impact == Impact.INFRAME || x.impact == Impact.SYNONYMOUS}
//        val dndsMutations = DndsMutation.fromFile(mutationsFile).filter { x -> x.impact == Impact.MISSENSE || x.impact == Impact.INFRAME }
//        val other = dndsMutations.filter { x -> x.impact != Impact.MISSENSE && x.impact != Impact.INFRAME }.sortedBy { x -> x.sample }


        logger.info("Loading dNdScv values: $dndsCVFile")
        val dndsCv = DndsCvGene.fromFile(dndsCVFile).associateBy { x -> x.gene }

        logger.info("Calculating gene mutation summary")
        val oncoGeneMutations = MutationsGene.oncoGeneMutations(dndsMutations).associateBy { x -> x.gene }
        val tsgGeneMutations = MutationsGene.tsgGeneMutations(dndsMutations).associateBy { x -> x.gene }

        // Log any missing genes from dNdS
        val dndsGenes = dndsCv.keys
        val missingDndsGenes = mutableSetOf<String>()
        missingDndsGenes.addAll(tsgGeneMutations.keys.subtract(dndsGenes))
        missingDndsGenes.addAll(oncoGeneMutations.keys.subtract(dndsGenes))
        for (gene in missingDndsGenes) {
            logger.warn("Dnds values missing for gene $gene")
        }

        val oncoLikelihood = LikelihoodGene(cohortSize, cohortLoad.totalLoad, dndsCv, oncoGeneMutations)
        val tsgLikelihood = LikelihoodGene(cohortSize, cohortLoad.totalLoad, dndsCv, tsgGeneMutations)
//        val tsgBiallelicLikelihood = LikelihoodGene(cohortSize, cohortLoad.biallelicLoad, dndsCv, tsgGeneMutations)
//        val tsgNonBiallelicLikelihood = LikelihoodGene(cohortSize, cohortLoad.nonBiallelicLoad, dndsCv, tsgGeneMutations)


        for (oncoGeneMutation in oncoLikelihood) {
            println(oncoGeneMutation)
        }

        for (tsgGeneMutation in tsgLikelihood) {
            println(tsgGeneMutation)
        }

    }

    private fun loadCohort(cohortFile: String): Pair<Int, TumorMutationalLoadSample> {
        val sampleLoads =  TumorMutationalLoadSample.fromFile(cohortFile)
        var biallelc = TumorMutationalLoad(0, 0, 0)
        var nonBiallelc = TumorMutationalLoad(0, 0, 0)
        for (sampleLoad in sampleLoads) {
            biallelc += sampleLoad.biallelicLoad
            nonBiallelc += sampleLoad.nonBiallelicLoad
        }

        return Pair(sampleLoads.size, TumorMutationalLoadSample("Total", biallelc, nonBiallelc, biallelc + nonBiallelc))
    }

    override fun close() {
        logger.info("Finished in ${(System.currentTimeMillis() - startTime) / 1000} seconds")
    }
}