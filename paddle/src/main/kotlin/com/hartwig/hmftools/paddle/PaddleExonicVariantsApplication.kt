package com.hartwig.hmftools.paddle

import com.hartwig.hmftools.common.drivercatalog.dnds.DndsVariantFile
import com.hartwig.hmftools.paddle.cohort.HighestPuritySample
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess
import org.apache.commons.cli.CommandLine
import org.apache.commons.cli.DefaultParser
import org.apache.commons.cli.Options
import org.apache.commons.cli.ParseException
import org.apache.logging.log4j.LogManager
import java.io.File

fun main(args: Array<String>) {

    fun createBasicOptions(): Options {
        val options = Options()
        DatabaseAccess.addDatabaseCmdLineArgs(options)
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
        val MAX_REPEAT_COUNT = 7
    }

    private val startTime = System.currentTimeMillis()
    private val dbAccess = DatabaseAccess.databaseAccess(cmd)

    override fun run() {
        val hpcFile = "/Users/jon/hmf/repos/hmftools/paddle/src/main/resources/highestPurityCohort.tsv"
        val exonicSomaticsFile = "/Users/jon/hmf/repos/hmftools/paddle/src/main/resources/HmfMutations.tsv"
        fun somaticFilename(x:String) = "/Users/jon/hmf/analysis/dnds/somatics/${x}.exonic.somatics.tsv"

        val highestPurityCohort = HighestPuritySample.readFile(hpcFile).take(10)
        for (sample in highestPurityCohort) {
            val somaticFilename =  somaticFilename(sample.sampleId)
            val somaticFile = File(somaticFilename)
            if (somaticFile.exists()) {
                logger.info("Skipping ${sample.sampleId}")
            } else {
                logger.info("Processing ${sample.sampleId}")
                val exonicSomatics = dbAccess.readDndsVariants(MAX_REPEAT_COUNT, sample.sampleId)
                DndsVariantFile.write(somaticFilename, exonicSomatics)

                val mutationalLoad = dbAccess.readDndsMutationLoad(sample.sampleId)
                println("Sdf")
            }
        }

        logger.info("Combining exonics")
        DndsVariantFile.writeHeader(exonicSomaticsFile)
        for (sample in highestPurityCohort) {
            val somaticFilename =  somaticFilename(sample.sampleId)
            val somaticFile = File(somaticFilename)
            if (!somaticFile.exists()) {
                throw IllegalStateException("Missing sample ${sample.sampleId}")
            }
            DndsVariantFile.append(exonicSomaticsFile, DndsVariantFile.read(somaticFilename))
        }

    }

    override fun close() {
        dbAccess.close()
        PaddleDndsApplication.logger.info("Finished in ${(System.currentTimeMillis() - startTime) / 1000} seconds")
    }
}

