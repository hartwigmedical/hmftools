package com.hartwig.hmftools.idgenerator

import com.hartwig.hmftools.extensions.cli.createCommandLine
import com.hartwig.hmftools.extensions.cli.createRunModeCommandLine
import com.hartwig.hmftools.idgenerator.extensions.readCurrentIds
import com.hartwig.hmftools.idgenerator.extensions.writeHmfIds
import org.apache.commons.cli.CommandLine
import org.apache.commons.cli.Option
import org.apache.commons.cli.OptionGroup
import org.apache.commons.cli.Options
import org.apache.logging.log4j.LogManager
import java.io.File

private val logger = LogManager.getLogger("HmfIdApplication")

fun main(args: Array<String>) {
    val cmd = createOptions().createRunModeCommandLine("hmf-id", args)
    when {
        cmd.hasOption(CREATE_SINGLE_HASH_MODE) -> {
            val singleHashCmd = singleHashModeOptions().createCommandLine("hmf-id", args)
            runHmfIdSingleHash(singleHashCmd)
        }
        cmd.hasOption(CREATE_IDS_MODE) -> {
            val createIdsCmd = createIdsModeOptions().createCommandLine("hmf-id", args)
            runHmfIdCreateIds(createIdsCmd)
        }
        cmd.hasOption(UPDATE_IDS_MODE) -> {
            val updateIdsCmd = updateIdsModeOptions().createCommandLine("hmf-id", args)
            runHmfIdUpdateIds(updateIdsCmd)
        }
        cmd.hasOption(ANONYMISE_IDS_MODE) -> {
            val anonymiseIdsCmd = anonymiseIdsModeOptions().createCommandLine("hmf-id", args)
            runHmfIdAnonymiseIds(anonymiseIdsCmd)
        }
    }
}

private fun createOptions(): Options {
    val options = Options()
    val inputModeOptionGroup = OptionGroup()
    inputModeOptionGroup.addOption(Option.builder(CREATE_SINGLE_HASH_MODE).required().desc("create single hash").build())
    inputModeOptionGroup.addOption(Option.builder(CREATE_IDS_MODE).required().desc("create hmf ids").build())
    inputModeOptionGroup.addOption(Option.builder(UPDATE_IDS_MODE).required().desc("update hmf ids").build())
    inputModeOptionGroup.addOption(Option.builder(ANONYMISE_IDS_MODE).required().desc("anonymise ids").build())
    inputModeOptionGroup.isRequired = true
    options.addOptionGroup(inputModeOptionGroup)
    return options
}

private fun singleHashModeOptions(): Options {
    val options = Options()
    options.addOption(Option.builder(CREATE_SINGLE_HASH_MODE).required().desc("create single hash").build())
    options.addOption(Option.builder(PASSWORD).required().hasArg().desc("password").build())
    options.addOption(Option.builder(PATIENT_ID).required().hasArg().desc("patient id to convert to hash").build())
    return options
}

private fun createIdsModeOptions(): Options {
    val options = Options()
    options.addOption(Option.builder(CREATE_IDS_MODE).required().desc("create hmf ids").build())
    options.addOption(Option.builder(PASSWORD).required().hasArg().desc("password").build())
    options.addOption(Option.builder(PATIENT_IDS_FILE).required().hasArg().desc("file containing a list of patients per line").build())
    options.addOption(Option.builder(OUTPUT_FILE).required().hasArg().desc("output file location").build())
    return options
}

private fun updateIdsModeOptions(): Options {
    val options = Options()
    options.addOption(Option.builder(UPDATE_IDS_MODE).required().desc("update hmf ids").build())
    options.addOption(Option.builder(PASSWORD).required().hasArg().desc("password").build())
    options.addOption(Option.builder(OLD_PASSWORD).optionalArg(true).hasArg().desc("password used to generate hashes in HMF ids file").build())
    options.addOption(Option.builder(PATIENT_IDS_FILE).required().hasArg().desc("file containing a list of patients per line").build())
    options.addOption(Option.builder(OUTPUT_FILE).required().hasArg().desc("output file location").build())
    return options
}

private fun anonymiseIdsModeOptions(): Options {
    val options = Options()
    options.addOption(Option.builder(ANONYMISE_IDS_MODE).required().desc("create hmf ids").build())
    options.addOption(Option.builder(PASSWORD).required().hasArg().desc("password").build())
    options.addOption(Option.builder(PATIENT_IDS_FILE).required().hasArg().desc("file containing a list of patients per line").build())
    return options
}

private fun runHmfIdSingleHash(cmd: CommandLine) {
    val generator = IdGenerator(cmd.getOptionValue(PASSWORD))
    val patientId = cmd.getOptionValue(PATIENT_ID)
    val hash = generator.hash(patientId)
    logger.info("Generated hash id for $patientId: $hash ")
}

private fun runHmfIdCreateIds(cmd: CommandLine) {
    val generator = IdGenerator(cmd.getOptionValue(PASSWORD))
    val patientIds = readPatientIdsFile(cmd.getOptionValue(PATIENT_IDS_FILE))
    val idMapping = generator.generateIds(patientIds)
    File(cmd.getOptionValue(OUTPUT_FILE)).writeHmfIds(idMapping.values)
    logger.info("Created hmf ids for ${idMapping.size} patients")
}

private fun runHmfIdUpdateIds(cmd: CommandLine) {
    val password = cmd.getOptionValue(PASSWORD)
    val generator = IdGenerator(password)
    val patientIds = readPatientIdsFile(cmd.getOptionValue(PATIENT_IDS_FILE))
    val currentIds = IdGenerator::class.java.getResource(CURRENT_IDS_FILE).openStream().readCurrentIds()
    val newIdMapping = generator.updateIds(cmd.getOptionValue(OLD_PASSWORD, password), patientIds, currentIds)
    File(cmd.getOptionValue(OUTPUT_FILE)).writeHmfIds(newIdMapping.values)
    logger.info("Updated hmf ids for ${newIdMapping.size} patients")
}

private fun runHmfIdAnonymiseIds(cmd: CommandLine) {
    val currentIds = IdGenerator::class.java.getResource(CURRENT_IDS_FILE).openStream().readCurrentIds()
    val generator = IdGenerator(cmd.getOptionValue(PASSWORD))
    val patientIds = readPatientIdsFile(cmd.getOptionValue(PATIENT_IDS_FILE))

    println("OriginalId,AnonymousId")
    patientIds.forEach { patient ->
        run {
            val hash = generator.hash(patient)
            val hmfId: HmfId? = currentIds.find { it.hash == hash }
            println(patient + "," + (hmfId?.patientId ?: "Unknown"))
        }
    }
}

private fun readPatientIdsFile(patientIdsFile: String): List<String> {
    return File(patientIdsFile).readLines().shuffled()
}