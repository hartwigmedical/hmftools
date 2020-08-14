package com.hartwig.hmftools.paddle

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

        val mutationsFile = "/Users/jon/hmf/repos/hmftools/paddle/src/main/resources/DndsMutations.AR.tsv"
        val dndsCVFile = "/Users/jon/hmf/repos/hmftools/paddle/src/main/resources/HmfRefCDSCv.tsv"

        logger.info("Loading mutations: $mutationsFile")
        val dndsMutations = DndsMutation.fromFile(mutationsFile)
//        val dndsMutations = DndsMutation.fromFile(mutationsFile).filter { x -> x.impact == Impact.MISSENSE || x.impact == Impact.INFRAME || x.impact == Impact.SYNONYMOUS}
//        val dndsMutations = DndsMutation.fromFile(mutationsFile).filter { x -> x.impact == Impact.MISSENSE || x.impact == Impact.INFRAME }
//        val other = dndsMutations.filter { x -> x.impact != Impact.MISSENSE && x.impact != Impact.INFRAME }.sortedBy { x -> x.sample }


        logger.info("Loading dNdScv values: $dndsCVFile")
        val dndsCv = DndsCvGene.fromFile(dndsCVFile).associateBy { x -> x.gene }

        logger.info("Calculating gene mutation summary")
        val oncoGeneMutations = MutationsGene.oncoGeneMutations(dndsMutations).associateBy { x -> x.gene }
        val tsgGeneMutations = MutationsGene.tsgGeneMutations(dndsMutations).associateBy { x -> x.gene }

        val dndsGenes = dndsCv.keys
        val oncoGenes = oncoGeneMutations.keys
        val tsgGenes = tsgGeneMutations.keys

        // Log any missing genes from dNdS
        val missingDndsGenes = mutableSetOf<String>()
        missingDndsGenes.addAll(tsgGenes.subtract(dndsGenes))
        missingDndsGenes.addAll(oncoGenes.subtract(dndsGenes))
        for (gene in missingDndsGenes) {
            logger.warn("Dnds values missing for gene $gene")
        }

        val oncoLikelihood = LikelihoodGene(dndsCv, oncoGeneMutations)
        val tsgLikelihood = LikelihoodGene(dndsCv, tsgGeneMutations)


        for (oncoGeneMutation in oncoLikelihood) {
            println(oncoGeneMutation)
        }

        for (tsgGeneMutation in tsgLikelihood) {
            println(tsgGeneMutation)
        }

    }

    override fun close() {
        logger.info("Finished in ${(System.currentTimeMillis() - startTime) / 1000} seconds")
    }
}