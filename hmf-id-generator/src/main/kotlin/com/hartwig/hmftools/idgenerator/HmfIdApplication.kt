package com.hartwig.hmftools.idgenerator

import com.hartwig.hmftools.common.amber.AmberAnonymous
import com.hartwig.hmftools.common.amber.AmberPatient
import com.hartwig.hmftools.extensions.cli.createCommandLine
import com.hartwig.hmftools.extensions.cli.options.HmfOptions
import com.hartwig.hmftools.extensions.cli.options.strings.InputOption
import com.hartwig.hmftools.extensions.cli.options.strings.RequiredInputOption
import com.hartwig.hmftools.extensions.cli.options.strings.RequiredOutputOption
import com.hartwig.hmftools.extensions.csv.CsvReader
import com.hartwig.hmftools.extensions.csv.CsvWriter
import com.hartwig.hmftools.idgenerator.anonymizedIds.HmfSampleIdCsv
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess
import org.apache.commons.cli.CommandLine
import org.apache.logging.log4j.LogManager
import kotlin.system.exitProcess

private val logger = LogManager.getLogger("HmfIdApplication")

fun main(args: Array<String>) {
    logger.info("Running id-generator using AMBER database tables")

    val hmfOptions = HmfOptions()
    DatabaseAccess.addDatabaseCmdLineArgs(hmfOptions.options, true)
    hmfOptions.add(RequiredInputOption(PASSWORD, "password"))
    hmfOptions.add(InputOption(NEW_PASSWORD, "password used to generate hashes in HMF ids file"))
    hmfOptions.add(RequiredInputOption(HASH_FILE_IN, "input hash file location"))
    hmfOptions.add(RequiredOutputOption(HASH_FILE_OUT, "output hash file location"))

    run(hmfOptions.createCommandLine("hmf-id", args))
}

private fun run(cmd: CommandLine) {
    val password = cmd.getOptionValue(PASSWORD)
    val newPassword = cmd.getOptionValue(NEW_PASSWORD, password)
    val databaseAccess = DatabaseAccess.databaseAccess(cmd)

    val amberPatients = databaseAccess.readAmberPatients()
    logger.info("Retrieved ${amberPatients.size} samples from amberPatient table")

    val rawMappings = databaseAccess.readAmberAnonymous().toSet()
    val existingMappings = rawMappings.filter { x -> !x.deleted() }.toSet()
    val deletedMappings = rawMappings.filter { x -> x.deleted() }.toSet()

    logger.info("Retrieved ${rawMappings.size} sample mappings from amberAnonymous table")
    if (rawMappings.isEmpty()) {
        logger.error("amberAnonymous table seems to be truncated. Exiting")
        exitProcess(1)
    }

    val deletedPatientInfo = deletedMappingPatientInfo(deletedMappings, amberPatients)
    if (deletedPatientInfo.isNotEmpty()) {
        for (sample in deletedPatientInfo) {
            logger.error("amberPatient table contains record for deleted sample $sample")
        }
        exitProcess(1)
    }

    val hashFileIn = cmd.getOptionValue(HASH_FILE_IN)
    val currentIds = CsvReader.readCSVByName<HmfSampleIdCsv>(hashFileIn).map { it.toHmfSampleId() }
    logger.info("Retrieved ${currentIds.size} sample hashes from ${hashFileIn}")

    logger.info("Processing samples")
    val amberAnonymizer = PatientAnonymizer(password, newPassword)
    val result = amberAnonymizer.anonymize(amberPatients, currentIds)
    val newMappings = AnonymizedRecord(newPassword, result, amberPatients.map { it.sample() }).map { x -> x.toAmberAnonymous() }

    val droppedMappings = existingMappings.subtract(newMappings)
    if (droppedMappings.isNotEmpty()) {
        for (missing in droppedMappings) {
            logger.error("Previous mapping ${missing} is no longer found")
        }
        exitProcess(1)
    }

    // Write to file and database
    val hashFileOut = cmd.getOptionValue(HASH_FILE_OUT)
    logger.info("Writing ${result.size} samples hashes to ${hashFileOut}")
    CsvWriter.writeCSV(result.map { it.toCsv() }, hashFileOut)

    val databaseOutput = newMappings + deletedMappings
    logger.info("Writing ${databaseOutput.size} sample mappings to database")
    databaseAccess.writeAmberAnonymous(databaseOutput)

    logger.info("Complete")
}

private fun deletedMappingPatientInfo(deletedMappings: Set<AmberAnonymous>, patients: List<AmberPatient>): Set<String> {
    val patientSet = patients.map { x -> x.sample() }.toSet()
    return deletedMappings.map { x -> x.sampleId() }.filter { x -> patientSet.contains(x) }.toSet()
}
