package com.hartwig.hmftools.paddle

import com.hartwig.hmftools.paddle.cohort.TumorMutationalLoad
import com.hartwig.hmftools.paddle.cohort.TumorMutationalLoadSample
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

    fun splitDnds(biallelicFile: String, nonBiallelicFile: String): Pair<Map<Gene, DndsCvGene>, Map<Gene, DndsCvGene>> {
        val dndsCvBiallelic = DndsCvGene.fromFileNoIndel(biallelicFile).associateBy { x -> x.gene }
        val dndsCvNonBiallelic = DndsCvGene.fromFileNoIndel(nonBiallelicFile).associateBy { x -> x.gene }

        val splitGenes = mutableSetOf<String>()
        for (gene in dndsCvBiallelic.keys) {
            val biallelic = dndsCvBiallelic[gene]!!
            val nonBiallelic = dndsCvNonBiallelic[gene]!!
            if (biallelic.missense.wCv > nonBiallelic.missense.wCv) {
                splitGenes.add(gene)
            }
        }

        return Pair(dndsCvBiallelic.filter {  it.key in splitGenes }, dndsCvNonBiallelic.filter {  it.key in splitGenes })
    }

    override fun run() {
        val cohortFile = "/Users/jon/hmf/repos/hmftools/paddle/src/main/resources/HmfTMB.tsv"
        val dndsCVFile = "/Users/jon/hmf/repos/hmftools/paddle/src/main/resources/HmfRefCDSCv.tsv"
        val dndsBiallelicCVFile = "/Users/jon/hmf/repos/hmftools/paddle/src/main/resources/HmfRefCDSCvBiallelic.tsv"
        val dndsNonBiallelicCVFile = "/Users/jon/hmf/repos/hmftools/paddle/src/main/resources/HmfRefCDSCvNonBiallelic.tsv"
        val mutationsFile = "/Users/jon/hmf/repos/hmftools/paddle/src/main/resources/DndsMutations.tsv"

        logger.info("Loading dNdScv values: $dndsCVFile")
        val dndsCv = DndsCvGene.fromFile(dndsCVFile).associateBy { x -> x.gene }
        val (dndsCvBiallelic, dndsCvNonBiallelic) = splitDnds(dndsBiallelicCVFile, dndsNonBiallelicCVFile)
        val tsSplitGenes = dndsCvBiallelic.keys

        logger.info("Loading cohort: $cohortFile")
        val (cohortSize, cohortLoad) = loadCohort(cohortFile)

        logger.info("Loading mutations: $mutationsFile")
        val dndsMutations = DndsMutation.fromFile(mutationsFile)
        val tsgBiallelicMutations = dndsMutations.filter { it.gene in tsSplitGenes && it.biallelic }
        val tsgNonBiallelicMutations = dndsMutations.filter { it.gene in tsSplitGenes && !it.biallelic  }

        logger.info("Calculating gene mutation summary")
        val oncoGeneMutations = MutationsGene.oncoGeneMutations(dndsMutations).associateBy { x -> x.gene }
        val tsgGeneMutations = MutationsGene.tsgGeneMutations(dndsMutations).associateBy { x -> x.gene }
        val tsgBiallelicGeneMutations = MutationsGene.tsgGeneMutations(tsgBiallelicMutations).associateBy { x -> x.gene }
        val tsgNonBiallelicGeneMutations = MutationsGene.tsgGeneMutations(tsgNonBiallelicMutations).associateBy { x -> x.gene }

//        // Log any missing genes from dNdS
//        val dndsGenes = dndsCv.keys
//        val missingDndsGenes = mutableSetOf<String>()
//        missingDndsGenes.addAll(tsgGeneMutations.keys.subtract(dndsGenes))
//        missingDndsGenes.addAll(oncoGeneMutations.keys.subtract(dndsGenes))
//        for (gene in missingDndsGenes) {
//            logger.warn("Dnds values missing for gene $gene")
//        }

        val oncoLikelihood = LikelihoodGene(cohortSize, cohortLoad.totalLoad, dndsCv, oncoGeneMutations)
        val tsgLikelihood = LikelihoodGene(cohortSize, cohortLoad.totalLoad, dndsCv, tsgGeneMutations)
        val tsgBiallelicLikelihood = LikelihoodGene(cohortSize, cohortLoad.biallelicLoad, dndsCvBiallelic, tsgBiallelicGeneMutations)
        val tsgNoneBiallelicLikelihood = LikelihoodGene(cohortSize, cohortLoad.nonBiallelicLoad, dndsCvNonBiallelic, tsgNonBiallelicGeneMutations)

        LikelihoodGene.writeFile(false, "/Users/jon/hmf/repos/hmftools/hmf-common/src/main/resources/dnds/DndsDriverLikelihoodOnco.tsv", oncoLikelihood.values)
        LikelihoodGene.writeFile(false, "/Users/jon/hmf/repos/hmftools/hmf-common/src/main/resources/dnds/DndsDriverLikelihoodTsg.tsv", tsgLikelihood.values)
        LikelihoodGene.writeFile(true, "/Users/jon/hmf/repos/hmftools/hmf-common/src/main/resources/dnds/DndsDriverLikelihoodTsgBiallelic.tsv", tsgBiallelicLikelihood.values)
        LikelihoodGene.writeFile(true, "/Users/jon/hmf/repos/hmftools/hmf-common/src/main/resources/dnds/DndsDriverLikelihoodTsgNonBiallelic.tsv", tsgNoneBiallelicLikelihood.values)

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