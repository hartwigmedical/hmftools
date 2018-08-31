package com.hartwig.hmftools.portal.converter

import com.google.common.annotations.VisibleForTesting
import com.google.common.collect.ArrayListMultimap
import com.google.common.collect.Multimap
import com.hartwig.hmftools.common.context.ProductionRunContextFactory
import com.hartwig.hmftools.common.context.RunContext
import com.hartwig.hmftools.extensions.cli.createCommandLine
import com.hartwig.hmftools.portal.converter.extensions.readSamplesData
import com.hartwig.hmftools.portal.converter.extensions.somaticVcfPath
import com.hartwig.hmftools.portal.converter.ids.extractPatientId
import com.hartwig.hmftools.portal.converter.records.ssm.SimpleSomaticMutation
import com.hartwig.hmftools.portal.converter.records.ssm.SimpleSomaticMutationMetadata
import org.apache.commons.cli.Option
import org.apache.commons.cli.Options
import org.apache.logging.log4j.LogManager
import java.io.File
import java.nio.file.Files
import java.nio.file.Paths
import java.util.stream.Collectors

private val logger = LogManager.getLogger("PortalDataConverterApplication")

fun main(args: Array<String>) {
    val options = createOptions()
    val cmd = options.createCommandLine("portal-data-converter", args)
    val patientDataFile = File(cmd.getOptionValue(CLINICAL_DATA_CSV))
    val runs = Files.lines(Paths.get(cmd.getOptionValue(RUNS_FILE))).map { it.trim() }.filter { it.isNotEmpty() }.collect(Collectors.toList()).map {
        ProductionRunContextFactory.fromRunDirectory(it)
    }
    convertSamples(runs, cmd.getOptionValue(OUTPUT_DIR), patientDataFile)
}

private fun createOptions(): Options {
    val options = Options()
    options.addOption(Option.builder(RUNS_FILE).required().hasArg().desc("path to file containing runs to convert").build())
    options.addOption(Option.builder(CLINICAL_DATA_CSV).required().hasArg().desc("path to csv file containing clinical data").build())
    options.addOption(Option.builder(OUTPUT_DIR).required().hasArg().desc("output directory").build())
    return options
}

@VisibleForTesting
internal fun convertSamples(runContexts: List<RunContext>, outputDirectory: String, patientDataFile: File) {
    val clinicalRecords: Multimap<String, SampleClinicalData> = ArrayListMultimap.create()
    val samplesData = patientDataFile.readSamplesData()
    val dataPerSample = samplesData.associateBy { it.sampleId }
    val dataPerPatient = samplesData.associateBy { it.cpctId }
    runContexts.forEach { it ->
        val (clinicalData, somaticVcfPath) = getPatientDataAndVcf(it, dataPerSample, dataPerPatient) ?: return@forEach
        if (clinicalData.tumorLocation == "Missing" || clinicalData.tumorLocation == "Other" || clinicalData.tumorLocation.toLowerCase().contains("unknown")) {
            logger.info("Skipping patient with ${clinicalData.tumorLocation} cancer type: ${clinicalData.hmfId}")
            return@forEach
        }
        val projectFolder = Files.createDirectories(Paths.get(outputDirectory + File.separator + "HMF-${clinicalData.tumorLocation}"))
        val folderPath = projectFolder.toAbsolutePath().toString()
        clinicalRecords.put(folderPath, clinicalData)
        TsvWriter.writeSomaticMutationMetadata(folderPath,
                                               clinicalData.sampleHmfId,
                                               listOf(SimpleSomaticMutationMetadata(clinicalData.sampleHmfId)))
        TsvWriter.writeSimpleSomaticMutation(folderPath,
                                             clinicalData.sampleHmfId,
                                             SimpleSomaticMutation(clinicalData.sampleHmfId, somaticVcfPath))
    }
    clinicalRecords.keySet().forEach { TsvWriter.writeClinicalData(it, clinicalRecords[it]) }
    logger.info("Done.")
}

private fun getPatientDataAndVcf(runContext: RunContext, dataPerSample: Map<String, SampleClinicalData>,
                                 dataPerPatient: Map<String, SampleClinicalData>): Pair<SampleClinicalData, String>? {
    val patientId = extractPatientId(runContext.tumorSample())
    val sampleClinicalData = dataPerSample[runContext.tumorSample()]
    val somaticVcfPath = runContext.somaticVcfPath()
    return when {
        patientId == null          -> {
            logger.warn("Could not extract valid patient id from sample: ${runContext.tumorSample()}")
            null
        }
        somaticVcfPath == null     -> {
            logger.warn("Could not locate somatic vcf file for set: ${runContext.setName()}")
            null
        }
        sampleClinicalData == null -> {
            val patientClinicalData = dataPerPatient[patientId]
            if (patientClinicalData == null) {
                logger.warn("Missing clinical data for patient: $patientId")
                null
            } else {
                Pair(patientClinicalData, somaticVcfPath)
            }
        }
        else                       -> Pair(sampleClinicalData, somaticVcfPath)
    }
}
