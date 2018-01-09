package com.hartwig.hmftools.portal.converter

import com.google.common.collect.ArrayListMultimap
import com.google.common.collect.Multimap
import com.hartwig.hmftools.common.context.ProductionRunContextFactory
import com.hartwig.hmftools.common.context.RunContext
import com.hartwig.hmftools.common.extensions.cli.createCommandLine
import com.hartwig.hmftools.portal.converter.extensions.readCpctSamples
import com.hartwig.hmftools.portal.converter.extensions.somaticVcfPath
import com.hartwig.hmftools.portal.converter.ids.extractPatientId
import com.hartwig.hmftools.portal.converter.records.SampleRecords
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
    val patientData = File(cmd.getOptionValue(CLINICAL_DATA_CSV)).readCpctSamples()
    val runs = Files.lines(Paths.get(cmd.getOptionValue(RUNS_FILE))).map { it.trim() }.filter { it.isNotEmpty() }
            .collect(Collectors.toList()).map { ProductionRunContextFactory.fromRunDirectory(it) }
    convertSamples(runs, cmd.getOptionValue(OUTPUT_DIR), patientData)
}

private fun createOptions(): Options {
    val options = Options()
    options.addOption(Option.builder(RUNS_FILE).required().hasArg().desc("path to file containing runs to convert").build())
    options.addOption(Option.builder(CLINICAL_DATA_CSV).required().hasArg()
            .desc("path to csv file containing clinical data").build())
    options.addOption(Option.builder(OUTPUT_DIR).required().hasArg().desc("output directory").build())
    return options
}

private fun convertSamples(runContexts: List<RunContext>, outputDirectory: String, patientData: Map<String, SampleClinicalData>) {
    val clinicalRecords: Multimap<String, SampleRecords> = ArrayListMultimap.create()

    runContexts.forEach { it ->
        val (patientId, clinicalData, somaticVcfPath) = getPatientDataAndVcf(it, patientData) ?: return
        val projectFolder = "$outputDirectory/HMF-${clinicalData.cancerType}"
        clinicalRecords.put(projectFolder, SampleRecords(patientId, clinicalData))
        TsvWriter.writeSomaticMutationMetadata(projectFolder, clinicalData.sampleId,
                listOf(SimpleSomaticMutationMetadata(clinicalData.sampleId)))
        TsvWriter.writeSimpleSomaticMutation(projectFolder, clinicalData.sampleId, SimpleSomaticMutation(somaticVcfPath))
    }
    clinicalRecords.keySet().forEach { TsvWriter.writeSampleRecords(it, clinicalRecords[it]) }
}

private fun getPatientDataAndVcf(runContext: RunContext,
                                 patientData: Map<String, SampleClinicalData>): Triple<String, SampleClinicalData, String>? {
    val clinicalData = patientData[runContext.tumorSample()]
    val somaticVcfPath = runContext.somaticVcfPath()
    val patientId = extractPatientId(runContext.tumorSample())
    return when {
        clinicalData == null -> {
            logger.warn("Missing clinical data for sample id: ${runContext.tumorSample()}")
            null
        }
        somaticVcfPath == null -> {
            logger.warn("Could not locate somatic vcf file for set: ${runContext.setName()}")
            null
        }
        patientId == null -> {
            logger.warn("Could not extract valid patient id from sample: ${runContext.tumorSample()}")
            null
        }
        else -> Triple(patientId, clinicalData, somaticVcfPath)
    }
}
