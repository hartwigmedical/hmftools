package com.hartwig.hmftools.paddle

import com.hartwig.hmftools.common.cli.Configs
import com.hartwig.hmftools.paddle.PaddleCohortApplication.Companion.MIN_PURITY
import com.hartwig.hmftools.paddle.PaddleCohortApplication.Companion.MIN_PURITY_DEFAULT
import com.hartwig.hmftools.paddle.PaddleCohortApplication.Companion.OUT
import com.hartwig.hmftools.paddle.PaddleCohortApplication.Companion.logger
import com.hartwig.hmftools.paddle.cohort.HighestPuritySample
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess
import org.apache.commons.cli.*
import org.apache.logging.log4j.LogManager

fun main(args: Array<String>) {
    fun createBasicOptions(): Options {
        val options = Options()
        val outOption = Option(OUT, true, "Path to output")
        outOption.isRequired = true
        options.addOption(outOption)

        val minPurityOption = Option(MIN_PURITY, true, "Min purity [$MIN_PURITY_DEFAULT]")
        options.addOption(minPurityOption)

        DatabaseAccess.addDatabaseCmdLineArgs(options)
        return options
    }

    @Throws(ParseException::class)
    fun createCommandLine(args: Array<String>, options: Options): CommandLine {
        return DefaultParser().parse(options, args)
    }
    val options = createBasicOptions()
    try {
        val cmd = createCommandLine(args, options)
        PaddleCohortApplication(cmd).use { x -> x.run() }
    } catch (e: ParseException) {
        logger.warn(e)
        val formatter = HelpFormatter()
        formatter.printHelp("PaddleCohort", options)
    } catch (e: Exception) {
        println(e)
    }
}

class PaddleCohortApplication(cmd: CommandLine) : AutoCloseable, Runnable {

    companion object {
        const val OUT = "out"
        const val MIN_PURITY = "min_purity"
        const val MIN_PURITY_DEFAULT = 0.2
        val logger = LogManager.getLogger(this::class.java)
    }

    private val minPurity = Configs.defaultDoubleValue(cmd, MIN_PURITY, MIN_PURITY_DEFAULT)
    private val outputFile = cmd.getOptionValue(OUT)
    private val startTime = System.currentTimeMillis()
    private val dbAccess = DatabaseAccess.databaseAccess(cmd)

    override fun run() {
        logger.info("Writing to $outputFile")
        val highestPurityCohort = HighestPuritySample.highestPurityCohort(minPurity, dbAccess)
        HighestPuritySample.writeFile(outputFile, highestPurityCohort)
    }

    override fun close() {
        dbAccess.close()
        PaddleDndsApplication.logger.info("Finished in ${(System.currentTimeMillis() - startTime) / 1000} seconds")
    }
}

