package com.hartwig.hmftools.idgenerator

import com.hartwig.hmftools.extensions.cli.createCommandLine
import com.hartwig.hmftools.extensions.cli.createRunModeCommandLine
import com.hartwig.hmftools.extensions.cli.options.HmfOptions
import com.hartwig.hmftools.extensions.cli.options.filesystem.RequiredInputFileOption
import com.hartwig.hmftools.extensions.cli.options.flags.RequiredFlagOption
import com.hartwig.hmftools.extensions.cli.options.strings.InputOption
import com.hartwig.hmftools.extensions.cli.options.strings.RequiredInputOption
import com.hartwig.hmftools.extensions.cli.options.strings.RequiredOutputOption
import com.hartwig.hmftools.extensions.csv.CsvReader
import com.hartwig.hmftools.extensions.csv.CsvWriter
import com.hartwig.hmftools.idgenerator.ids.PatientId
import com.hartwig.hmftools.idgenerator.ids.SampleId
import org.apache.commons.cli.CommandLine
import org.apache.commons.cli.Option
import org.apache.commons.cli.OptionGroup
import org.apache.commons.cli.Options
import org.apache.logging.log4j.LogManager
import java.io.File
import java.util.*

private val logger = LogManager.getLogger("HmfIdApplication")

fun main(args: Array<String>) {
    logger.info("Running id-generator $Version")
    val cmd = createOptions().createRunModeCommandLine("hmf-id", args)
    when {
        cmd.hasOption(CREATE_SINGLE_HASH_MODE) -> {
            val singleHashCmd = singleHashModeOptions().createCommandLine("hmf-id", args)
            runSingleHash(singleHashCmd)
        }
        cmd.hasOption(CREATE_IDS_MODE)         -> {
            val createIdsCmd = createIdsModeOptions().createCommandLine("hmf-id", args)
            runCreateIds(createIdsCmd)
        }
        cmd.hasOption(UPDATE_IDS_MODE)         -> {
            val updateIdsCmd = updateIdsModeOptions().createCommandLine("hmf-id", args)
            runUpdateIds(updateIdsCmd)
        }
        cmd.hasOption(ANONYMIZE_IDS_MODE)      -> {
            val anonymizeIdsCmd = anonymizeIdsModeOptions().createCommandLine("hmf-id", args)
            runAnonymizeIds(anonymizeIdsCmd)
        }
    }
}

private fun createOptions(): Options {
    val options = Options()
    val inputModeOptionGroup = OptionGroup()
    inputModeOptionGroup.addOption(Option.builder(CREATE_SINGLE_HASH_MODE).required().desc("create single hash").build())
    inputModeOptionGroup.addOption(Option.builder(CREATE_IDS_MODE).required().desc("create hmf ids").build())
    inputModeOptionGroup.addOption(Option.builder(UPDATE_IDS_MODE).required().desc("update hmf ids").build())
    inputModeOptionGroup.addOption(Option.builder(ANONYMIZE_IDS_MODE).required().desc("anonymize ids").build())
    inputModeOptionGroup.isRequired = true
    options.addOptionGroup(inputModeOptionGroup)
    return options
}

private fun singleHashModeOptions(): HmfOptions {
    val hmfOptions = HmfOptions()
    hmfOptions.add(RequiredFlagOption(CREATE_SINGLE_HASH_MODE, "create single hash"))
    hmfOptions.add(RequiredInputOption(PASSWORD, "password"))
    hmfOptions.add(RequiredInputOption(SAMPLE_ID, "sample id to convert to hash"))
    return hmfOptions
}

private fun createIdsModeOptions(): HmfOptions {
    val hmfOptions = HmfOptions()
    hmfOptions.add(RequiredFlagOption(CREATE_IDS_MODE, "create hmf ids"))
    hmfOptions.add(RequiredInputOption(PASSWORD, "password"))
    hmfOptions.add(RequiredInputFileOption(SAMPLE_IDS_FILE, "file containing a list of samples, one per line"))
    hmfOptions.add(RequiredInputFileOption(PATIENT_MAPPING_FILE, "csv containing the patient mapping, a patient pair per line"))
    hmfOptions.add(RequiredOutputOption(OUTPUT_FILE, "output file location"))
    hmfOptions.add(RequiredOutputOption(SAMPLE_MAPPING_OUTPUT_FILE, "sample mapping output file location"))
    return hmfOptions
}

private fun updateIdsModeOptions(): HmfOptions {
    val hmfOptions = HmfOptions()
    hmfOptions.add(RequiredFlagOption(UPDATE_IDS_MODE, "update hmf ids"))
    hmfOptions.add(RequiredInputOption(PASSWORD, "password"))
    hmfOptions.add(InputOption(NEW_PASSWORD, "password used to generate hashes in HMF ids file"))
    hmfOptions.add(RequiredInputFileOption(SAMPLE_IDS_FILE, "file containing a list of samples, one per line"))
    hmfOptions.add(RequiredInputFileOption(PATIENT_MAPPING_FILE, "csv containing the patient mapping, a patient pair per line"))
    hmfOptions.add(RequiredOutputOption(OUTPUT_FILE, "output file location"))
    hmfOptions.add(RequiredOutputOption(SAMPLE_MAPPING_OUTPUT_FILE, "sample mapping output file location"))
    return hmfOptions
}

private fun anonymizeIdsModeOptions(): HmfOptions {
    val hmfOptions = HmfOptions()
    hmfOptions.add(RequiredFlagOption(ANONYMIZE_IDS_MODE, "anonymize ids"))
    hmfOptions.add(RequiredInputOption(PASSWORD, "password"))
    hmfOptions.add(RequiredInputFileOption(SAMPLE_IDS_FILE, "file containing a list of samples, one per line"))
    hmfOptions.add(RequiredInputFileOption(PATIENT_MAPPING_FILE, "csv containing the patient mapping, a patient pair per line"))
    return hmfOptions
}

private fun runSingleHash(cmd: CommandLine) {
    logger.info("Mode: single hash")
    val generator = IdGenerator(cmd.getOptionValue(PASSWORD))
    val sampleId = SampleId(cmd.getOptionValue(SAMPLE_ID))
    sampleId ?: return
    val hash = generator.hash(sampleId.id)
    logger.info("Generated hash id for ${sampleId.id}: ${hash.value} ")
}

private fun runCreateIds(cmd: CommandLine) {
    logger.info("Mode: create ids")
    val password = cmd.getOptionValue(PASSWORD)
    val anonymizer = SampleAnonymizer(password)
    val samplesInput = readSamplesInput(cmd.getOptionValue(SAMPLE_IDS_FILE), cmd.getOptionValue(PATIENT_MAPPING_FILE))
    val anonymizedSamples = anonymizer.anonymize(password, samplesInput, emptyList())
    val sampleMappingRecords = anonymizedSamples.sampleMapping.map { HmfSampleMappingRecord(it.key, it.value) }
    CsvWriter.writeCSV(anonymizedSamples.map { HmfSampleIdRecord(it) }, cmd.getOptionValue(OUTPUT_FILE))
    CsvWriter.writeCSV(sampleMappingRecords, cmd.getOptionValue(SAMPLE_MAPPING_OUTPUT_FILE))
    logger.info("Created hmf sample ids for ${anonymizedSamples.size} samples")
}

private fun runUpdateIds(cmd: CommandLine) {
    logger.info("Mode: update ids")
    val password = cmd.getOptionValue(PASSWORD)
    val newPassword = cmd.getOptionValue(NEW_PASSWORD, password)
    val anonymizer = SampleAnonymizer(password)
    val samplesInput = readSamplesInput(cmd.getOptionValue(SAMPLE_IDS_FILE), cmd.getOptionValue(PATIENT_MAPPING_FILE))
    val currentIds = CsvReader.readCSVByName<HmfSampleIdRecord>(IdGenerator::class.java.getResource(CURRENT_IDS_FILE).openStream())
            .map { it.toHmfSampleId() }
    logger.info("Read ${currentIds.size} anonymized samples")

    val anonymizedSamples = anonymizer.anonymize(newPassword, samplesInput, currentIds)
    val sampleMappingRecords = anonymizedSamples.sampleMapping.map { HmfSampleMappingRecord(it.key, it.value) }
    CsvWriter.writeCSV(anonymizedSamples.map { HmfSampleIdRecord(it) }, cmd.getOptionValue(OUTPUT_FILE))
    CsvWriter.writeCSV(sampleMappingRecords, cmd.getOptionValue(SAMPLE_MAPPING_OUTPUT_FILE))
    logger.info("Created hmf sample ids for ${anonymizedSamples.size} samples")
}

private fun runAnonymizeIds(cmd: CommandLine) {
    logger.info("Mode: anonymize ids")
    val currentIds = CsvReader.readCSVByName<HmfSampleIdRecord>(IdGenerator::class.java.getResource(CURRENT_IDS_FILE).openStream())
            .map { it.toHmfSampleId() }
    val samplesInput = readSamplesInput(cmd.getOptionValue(SAMPLE_IDS_FILE), cmd.getOptionValue(PATIENT_MAPPING_FILE))
    val anonymizedSamples = AnonymizedSamples(cmd.getOptionValue(PASSWORD), currentIds, samplesInput)
    println("OriginalId,AnonymousId")
    samplesInput.samples.sortedWith(Comparator.comparing<SampleId, String> { it.id }).forEach { sample ->
        anonymizedSamples[sample]
        println("${sample.id},${anonymizedSamples[sample]?.plaintext ?: "Unknown"}")
    }
}

private fun readSamplesInput(sampleIdsFile: String, patientMappingFile: String): SamplesInput {
    val samples = File(sampleIdsFile).readLines().shuffled().mapNotNull { SampleId(it) }
    val patientMapping = File(patientMappingFile).readLines().mapNotNull {
        val ids = it.split(",")
        val patientId = PatientId(ids[0])
        val canonicalId = PatientId(ids[1])
        return@mapNotNull if (patientId == null || canonicalId == null) null else patientId to canonicalId
    }.toMap()
    return SamplesInput(samples, patientMapping)
}
