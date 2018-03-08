package com.hartwig.hmftools.idgenerator

import com.hartwig.hmftools.common.ecrf.projections.ImmutablePortalClinicalData
import com.hartwig.hmftools.common.ecrf.projections.PortalClinicalData
import com.hartwig.hmftools.extensions.cli.createCommandLine
import com.hartwig.hmftools.extensions.cli.createRunModeCommandLine
import com.hartwig.hmftools.idgenerator.extensions.readOldIds
import com.hartwig.hmftools.idgenerator.extensions.writeHmfIds
import org.apache.commons.cli.CommandLine
import org.apache.commons.cli.Option
import org.apache.commons.cli.OptionGroup
import org.apache.commons.cli.Options
import org.apache.logging.log4j.LogManager
import java.io.File
import java.io.IOException
import java.nio.file.Files
import java.nio.file.Paths

private val logger = LogManager.getLogger("HmfIdApplication")

fun main(args: Array<String>) {
    val cmd = createOptions().createRunModeCommandLine("hmf-id", args)
    when {
        cmd.hasOption(CREATE_MODE) -> {
            val createCmd = createModeOptions().createCommandLine("hmf-id", args)
            runHmfIdCreate(createCmd)
        }
        cmd.hasOption(UPDATE_MODE) -> {
            val updateCmd = updateModeOptions().createCommandLine("hmf-id", args)
            runHmfIdUpdate(updateCmd)
        }
    }
}

private fun createOptions(): Options {
    val options = Options()
    val inputModeOptionGroup = OptionGroup()
    inputModeOptionGroup.addOption(Option.builder(CREATE_MODE).required().desc("create hmf ids").build())
    inputModeOptionGroup.addOption(Option.builder(UPDATE_MODE).required().desc("update hmf ids").build())
    inputModeOptionGroup.isRequired = true
    options.addOptionGroup(inputModeOptionGroup)
    return options
}

private fun createModeOptions(): Options {
    val options = Options()
    options.addOption(Option.builder(CREATE_MODE).required().desc("create hmf ids").build())
    options.addOption(Option.builder(PASSWORD).required().hasArg().desc("password").build())
    options.addOption(Option.builder(PORTAL_CLINICAL_DATA).required().hasArg().desc("CSV file containing portal clinical data dump").build())
    options.addOption(Option.builder(OUTPUT_FILE).required().hasArg().desc("output file location").build())
    return options
}

private fun updateModeOptions(): Options {
    val options = Options()
    options.addOption(Option.builder(UPDATE_MODE).required().desc("update hmf ids").build())
    options.addOption(Option.builder(PASSWORD).required().hasArg().desc("password").build())
    options.addOption(Option.builder(OLD_PASSWORD).required().hasArg().desc("password used to generate hashes in HMF ids file").build())
    options.addOption(Option.builder(PORTAL_CLINICAL_DATA).required().hasArg().desc("CSV file containing portal clinical data dump").build())
    options.addOption(Option.builder(HMF_IDS_FILE).required().hasArg().desc("file containing current mapping of patients to HMF ids").build())
    options.addOption(Option.builder(OUTPUT_FILE).required().hasArg().desc("output file location").build())
    return options
}

private fun runHmfIdCreate(cmd: CommandLine) {
    val generator = IdGenerator(cmd.getOptionValue(PASSWORD))
    val clinicalData = PortalClinicalData.readRecords(cmd.getOptionValue(PORTAL_CLINICAL_DATA))
    val idMapping = generator.generateIds(clinicalData.map { it.cpctId() })
    writeAnonymizedPatients(File(cmd.getOptionValue(PORTAL_CLINICAL_DATA)), clinicalData, idMapping)
    File(cmd.getOptionValue(OUTPUT_FILE)).writeHmfIds(idMapping.values)
}

private fun runHmfIdUpdate(cmd: CommandLine) {
    val generator = IdGenerator(cmd.getOptionValue(PASSWORD))
    val clinicalData = PortalClinicalData.readRecords(cmd.getOptionValue(PORTAL_CLINICAL_DATA))
    val oldIds = File(cmd.getOptionValue(HMF_IDS_FILE)).readOldIds()
    val newIdMapping = generator.updateIds(cmd.getOptionValue(OLD_PASSWORD), clinicalData.map { it.cpctId() }, oldIds)
    writeAnonymizedPatients(File(cmd.getOptionValue(PORTAL_CLINICAL_DATA)), clinicalData, newIdMapping)
    File(cmd.getOptionValue(OUTPUT_FILE)).writeHmfIds(newIdMapping.values)
}

private fun writeAnonymizedPatients(patientsFile: File, clinicalData: List<PortalClinicalData>, hmfIdMapping: Map<CpctId, HmfId>) {
    val canonicalPatientsFile = patientsFile.canonicalFile
    val anonymizedData = clinicalData.map { anonymizePatient(it, hmfIdMapping) }
    val outputPath = if (Files.isSymbolicLink(canonicalPatientsFile.toPath())) {
        val symlinkTarget = Files.readSymbolicLink(canonicalPatientsFile.toPath())
        symlinkTarget.toFile().parent + File.separator + "anonymized_${symlinkTarget.fileName.toFile().name}"
    } else {
        canonicalPatientsFile.parent + File.separator + "anonymized_${canonicalPatientsFile.name}"
    }
    PortalClinicalData.writeRecords(outputPath, anonymizedData)
    updateAnonymizedSymlink(canonicalPatientsFile, outputPath)
}

private fun updateAnonymizedSymlink(patientsFile: File, outputPath: String) {
    if (Files.isSymbolicLink(patientsFile.toPath())) {
        val linkPath = Paths.get(patientsFile.parent + File.separator + "anonymized_${patientsFile.name}")
        try {
            Files.deleteIfExists(linkPath)
            Files.createSymbolicLink(linkPath, Paths.get(outputPath))
        } catch (e: IOException) {
            logger.warn("Failed to update symlink {}. Cause: {}", linkPath, e.message)
        }
    }
}

private fun anonymizePatient(patientData: PortalClinicalData, hmfIdMapping: Map<CpctId, HmfId>): PortalClinicalData {
    val patientHmfId = "HMF${hmfIdMapping[patientData.cpctId()]!!.id}"
    val sampleHmfId = if (patientData.sampleId().length > 12 && patientData.sampleId()[12] == 'T') {
        patientHmfId + patientData.sampleId().substring(12)
    } else {
        patientData.sampleId()
    }
    return ImmutablePortalClinicalData.builder().from(patientData).hmfId(patientHmfId).sampleHmfId(sampleHmfId).build()
}
