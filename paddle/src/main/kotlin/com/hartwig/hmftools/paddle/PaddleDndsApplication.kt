package com.hartwig.hmftools.paddle

import com.hartwig.hmftools.paddle.PaddleDndsApplication.Companion.WORK_DIR
import com.hartwig.hmftools.paddle.cohort.CohortLoad
import com.hartwig.hmftools.paddle.dnds.DndsCvGene
import com.hartwig.hmftools.paddle.dnds.DndsMutation
import com.hartwig.hmftools.paddle.likelihood.LikelihoodGene
import com.hartwig.hmftools.paddle.mutation.MutationsGene
import org.apache.commons.cli.*
import org.apache.logging.log4j.LogManager

fun main(args: Array<String>) {

    fun createBasicOptions(): Options {
        val options = Options()

        val outputDirOption = Option(WORK_DIR, true, "Directory where to find all the inputs fro Dnds")
        outputDirOption.isRequired = true
        options.addOption(outputDirOption)

        return options
    }

    @Throws(ParseException::class)
    fun createCommandLine(args: Array<String>, options: Options): CommandLine {
        return DefaultParser().parse(options, args)
    }

    try {
        val options = createBasicOptions()
        val cmd = createCommandLine(args, options)

        PaddleDndsApplication(cmd).use { x -> x.run() }
    } catch (e: Exception) {
        println(e)
    }
}

class PaddleDndsApplication(cmd: CommandLine) : AutoCloseable, Runnable {

    companion object {
        val logger = LogManager.getLogger(this::class.java)

        const val WORK_DIR = "work_dir"
    }

    private val startTime = System.currentTimeMillis()
    private val workDir = cmd.getOptionValue(WORK_DIR)

    override fun run() {
        val cohortFile = "${workDir}/mutationalLoad.tsv"
        val dndsCVFile = "${workDir}/HmfRefCDSCv.tsv"
        val mutationsFile = "${workDir}/DndsMutations.tsv"

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

        LikelihoodGene.writeFile(false, "${workDir}/DndsDriverLikelihoodOnco.tsv", oncoLikelihood.values)
        LikelihoodGene.writeFile(false, "${workDir}/DndsDriverLikelihoodTsg.tsv", tsgLikelihood.values)
    }

    override fun close() {
        logger.info("Finished in ${(System.currentTimeMillis() - startTime) / 1000} seconds")
    }
}