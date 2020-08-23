package com.hartwig.hmftools.paddle

import com.hartwig.hmftools.paddle.cohort.HighestPuritySample
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess
import org.apache.commons.cli.CommandLine
import org.apache.commons.cli.DefaultParser
import org.apache.commons.cli.Options
import org.apache.commons.cli.ParseException
import org.apache.logging.log4j.LogManager

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

        PaddleCohortApplication(cmd).use { x -> x.run() }
    } catch (e: Exception) {
        println(e)
    }
}

class PaddleCohortApplication(cmd: CommandLine) : AutoCloseable, Runnable {

    companion object {
        val logger = LogManager.getLogger(this::class.java)
    }

    private val startTime = System.currentTimeMillis()
    private val dbAccess = DatabaseAccess.databaseAccess(cmd)


    override fun run() {
        val highestPurityCohort = HighestPuritySample.highestPurityCohort(0.2, dbAccess)
        HighestPuritySample.writeFile("/Users/jon/hmf/repos/hmftools/paddle/src/main/resources/highestPurityCohort.tsv", highestPurityCohort)
    }

    override fun close() {
        dbAccess.close()
        PaddleDndsApplication.logger.info("Finished in ${(System.currentTimeMillis() - startTime) / 1000} seconds")
    }
}

