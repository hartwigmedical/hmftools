package com.hartwig.hmftools.paddle

import com.hartwig.hmftools.common.drivercatalog.dnds.DndsMutationalLoadFile
import com.hartwig.hmftools.common.drivercatalog.dnds.DndsVariantFile
import com.hartwig.hmftools.paddle.PaddleExonicVariantsApplication.Companion.COHORT_TSV
import com.hartwig.hmftools.paddle.PaddleExonicVariantsApplication.Companion.OUTPUT_DIR
import com.hartwig.hmftools.paddle.cohort.CohortSample
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess
import org.apache.commons.cli.*
import org.apache.logging.log4j.LogManager
import java.io.File

fun main(args: Array<String>) {

    fun createBasicOptions(): Options {
        val options = Options()
        DatabaseAccess.addDatabaseCmdLineArgs(options)

        val outputDirOption = Option(OUTPUT_DIR, true, "Directory in which to write the output")
        outputDirOption.isRequired = true
        options.addOption(outputDirOption)

        val cohortTsvOption = Option(COHORT_TSV, true, "TSV file which contains the cohort to run on")
        cohortTsvOption.isRequired = true
        options.addOption(cohortTsvOption)

        return options
    }

    @Throws(ParseException::class)
    fun createCommandLine(args: Array<String>, options: Options): CommandLine {
        return DefaultParser().parse(options, args)
    }

    try {
        val options = createBasicOptions()
        val cmd = createCommandLine(args, options)

        PaddleExonicVariantsApplication(cmd).use { x -> x.run() }
    } catch (e: Exception) {
        println(e)
    }
}

class PaddleExonicVariantsApplication(cmd: CommandLine) : AutoCloseable, Runnable {

    companion object {
        val logger = LogManager.getLogger(this::class.java)
        const val MAX_REPEAT_COUNT = 7
        const val OUTPUT_DIR = "output_dir"
        const val COHORT_TSV = "cohort_tsv"
    }

    private val startTime = System.currentTimeMillis()
    private val dbAccess = DatabaseAccess.databaseAccess(cmd)
    private val outputDir = cmd.getOptionValue(OUTPUT_DIR)
    private val cohortFile = cmd.getOptionValue(COHORT_TSV)

    override fun run() {
        val cohortMutationalLoadFile = "$outputDir/mutationalLoad.tsv"
        val somaticsDir = "$outputDir/somatics"

        File(somaticsDir).mkdirs()

        fun somaticFilename(x:String) = "${somaticsDir}/${x}.exonic.somatics.tsv"

        if (!File(cohortMutationalLoadFile).exists()) {
            DndsMutationalLoadFile.writeHeader(cohortMutationalLoadFile)
        }

        val oldCohortMutationalLoad = DndsMutationalLoadFile.read(cohortMutationalLoadFile).associateBy { x -> x.sampleId() }

        val cohort = CohortSample.readFile(cohortFile)
        for (sample in cohort) {
            val somaticFilename =  somaticFilename(sample.sampleId)
            val somaticFile = File(somaticFilename)
            if (!somaticFile.exists()) {
                logger.info("Processing ${sample.sampleId}")
                val exonicSomatics = dbAccess.readDndsVariants(MAX_REPEAT_COUNT, sample.sampleId)
                DndsVariantFile.write(somaticFilename, exonicSomatics)
            }

            if (!oldCohortMutationalLoad.containsKey(sample.sampleId)) {
                val sampleMutationalLoad = dbAccess.readDndsMutationLoad(sample.sampleId)
                DndsMutationalLoadFile.append(cohortMutationalLoadFile, listOf(sampleMutationalLoad))
            }
        }
    }

    override fun close() {
        dbAccess.close()
        logger.info("Finished in ${(System.currentTimeMillis() - startTime) / 1000} seconds")
    }
}

