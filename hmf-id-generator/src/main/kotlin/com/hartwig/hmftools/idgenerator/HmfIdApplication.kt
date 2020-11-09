package com.hartwig.hmftools.idgenerator

import com.hartwig.hmftools.common.amber.AmberPatient
import com.hartwig.hmftools.extensions.cli.createCommandLine
import com.hartwig.hmftools.extensions.cli.options.HmfOptions
import com.hartwig.hmftools.extensions.cli.options.strings.InputOption
import com.hartwig.hmftools.extensions.cli.options.strings.RequiredInputOption
import com.hartwig.hmftools.extensions.cli.options.strings.RequiredOutputOption
import com.hartwig.hmftools.extensions.csv.CsvReader
import com.hartwig.hmftools.extensions.csv.CsvWriter
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess
import org.apache.commons.cli.CommandLine
import org.apache.logging.log4j.LogManager
import org.apache.logging.log4j.util.Strings
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
    val oldPassword = cmd.getOptionValue(PASSWORD)
    val newPassword = cmd.getOptionValue(NEW_PASSWORD, oldPassword)
    val databaseAccess = DatabaseAccess.databaseAccess(cmd)
    val hashFileIn = cmd.getOptionValue(HASH_FILE_IN)
    val oldGenerator = IdGenerator(oldPassword)
    val newGenerator = IdGenerator(newPassword)

    // Retrieve Data
    val amberPatients = databaseAccess.readAmberPatients()
    logger.info("Retrieved ${amberPatients.size} samples from amberPatient table")

    val currentDatabaseAnonymous = databaseAccess.readAmberAnonymous().map { x -> HmfSample(x) }
    logger.info("Retrieved ${currentDatabaseAnonymous.size} sample mappings from amberAnonymous table")

    val currentFileAnonymous = CsvReader.readCSVByName<HmfSampleCsv>(hashFileIn).toSet()
    logger.info("Retrieved ${currentFileAnonymous.size} sample mappings from $hashFileIn")

    // Validate amberPatient and amberAnonymous tables
    if (!validAmberPatients(currentDatabaseAnonymous, amberPatients)) {
        logger.error("amberPatient and amberAnonymous are not synchronized")
        exitProcess(1)
    }

    // Validate synchronization between database and file
    val currentDatabaseCsv = currentDatabaseAnonymous.map { x -> x.toCsv(oldGenerator) }.toSet()
    if (!databaseAndFileInSync(currentDatabaseCsv, currentFileAnonymous)) {
        logger.error("Database and file are not synchronized")
        exitProcess(1)
    }

    logger.info("Processing samples")
    val amberAnonymizer = PatientAnonymizer()
    val result = amberAnonymizer.anonymize(amberPatients, currentDatabaseAnonymous.toList())


    // Write to file and database
    val hashFileOut = cmd.getOptionValue(HASH_FILE_OUT)
    logger.info("Writing ${result.size} samples hashes to ${hashFileOut}")
    CsvWriter.writeCSV(result.map { x -> x.toCsv(newGenerator) }, hashFileOut)

    logger.info("Writing ${result.size} sample mappings to database")
    databaseAccess.writeAmberAnonymous(result.map { x -> x.toAmberAnonymous() })

    logger.info("Complete")
}


private fun validAmberPatients(currentDatabase: List<HmfSample>, patients: List<AmberPatient>): Boolean {
    val actualSamples = patients.map { x -> x.sample() }.toSet()
    val expectedSamples = currentDatabase.filter { x -> !x.deleted }.map { x -> x.sample }
    val missingSamples = expectedSamples subtract actualSamples
    for (missingSample in missingSamples) {
        logger.error("Missing sample $missingSample from amberPatient table")
    }

    val deletedSamples = currentDatabase.filter { x -> x.deleted }.map { x -> x.sample }
    val unexpectedSamples = actualSamples intersect deletedSamples
    for (unexpectedSample in unexpectedSamples) {
        logger.error("Deleted sample $unexpectedSample in amberPatient table")
    }

    if (missingSamples.isNotEmpty() || unexpectedSamples.isNotEmpty()) {
        return false
    }

    return true
}


private fun databaseAndFileInSync(currentDatabaseCsv: Set<HmfSampleCsv>, currentFileCsv: Set<HmfSampleCsv>): Boolean {
    // Ignore deleted flag for moment
    val adjustedCurrentDatabaseCsv = currentDatabaseCsv.map { x -> x.copy(deleted = false.toString()) }
    val adjustedCurrentFileCsv = currentFileCsv.map { x -> x.copy(deleted = false.toString()) }

    val missingFromDatabase = adjustedCurrentFileCsv subtract adjustedCurrentDatabaseCsv
    val missingFromFile = adjustedCurrentDatabaseCsv subtract adjustedCurrentFileCsv
    if (missingFromDatabase.isNotEmpty()) {
        for (missing in missingFromDatabase) {
            logger.error("${missing.hmfSampleId} missing from database")
        }
    }

    if (missingFromFile.isNotEmpty()) {
        for (missing in missingFromFile) {
            logger.error("${missing.hmfSampleId} missing from file")
        }
    }

    if (missingFromDatabase.isNotEmpty() || missingFromFile.isNotEmpty()) {
        return false
    }

    return true
}

private fun sample(password: String, hash: String): String {
    val generator = IdGenerator(password)
    val prefixes = setOf("WIDE", "CPCT", "DRUP")
    val locations = setOf("01", "02")
    val suffixes: Set<String> = setOf("T", "TI", "TII", "TIII")

    logger.info("Attempting to crack sample")

    for (prefix in prefixes) {
        logger.info("  using prefix $prefix")
        for (location in locations) {
            logger.info("   using location $location")
            for (suffix in suffixes) {
                logger.info("     using suffix $suffix")
                for (i in 1..500000) {
                    val sample = prefix + location + i.toString().padStart(6, '0') + suffix
                    val newHash = generator.hash(sample)
                    if (newHash == hash) {
                        logger.info("Success: $sample")
                        return sample
                    }
                }
            }

        }
    }
    return Strings.EMPTY;
}

