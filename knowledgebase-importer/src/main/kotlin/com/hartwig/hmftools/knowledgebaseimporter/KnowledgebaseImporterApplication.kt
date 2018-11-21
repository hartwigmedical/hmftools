package com.hartwig.hmftools.knowledgebaseimporter

import com.hartwig.hmftools.apiclients.iclusion.api.IclusionApiWrapper
import com.hartwig.hmftools.apiclients.iclusion.data.IclusionStudyDetails
import com.hartwig.hmftools.extensions.cli.createCommandLine
import com.hartwig.hmftools.extensions.cli.options.HmfOptions
import com.hartwig.hmftools.extensions.cli.options.commands.RequiredCommandOption
import com.hartwig.hmftools.extensions.cli.options.filesystem.RequiredInputFileOption
import com.hartwig.hmftools.extensions.cli.options.strings.RequiredInputOption
import com.hartwig.hmftools.extensions.cli.options.strings.RequiredOutputOption
import com.hartwig.hmftools.extensions.csv.CsvReader
import com.hartwig.hmftools.extensions.csv.CsvWriter
import com.hartwig.hmftools.knowledgebaseimporter.cgi.Cgi
import com.hartwig.hmftools.knowledgebaseimporter.civic.Civic
import com.hartwig.hmftools.knowledgebaseimporter.cosmic.Cosmic
import com.hartwig.hmftools.knowledgebaseimporter.diseaseOntology.DiseaseOntology
import com.hartwig.hmftools.knowledgebaseimporter.diseaseOntology.Doid
import com.hartwig.hmftools.knowledgebaseimporter.iclusion.Iclusion
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.RecordAnalyzer
import com.hartwig.hmftools.knowledgebaseimporter.oncoKb.OncoKb
import com.hartwig.hmftools.knowledgebaseimporter.output.CancerTypeDoidOutput
import com.hartwig.hmftools.knowledgebaseimporter.output.HmfDrug
import htsjdk.samtools.reference.IndexedFastaSequenceFile
import org.apache.commons.cli.CommandLine
import org.apache.logging.log4j.LogManager
import java.io.File

private val logger = LogManager.getLogger("KnowledgebaseImporterApplication")

fun main(args: Array<String>) {
    logger.info("Running Knowledgebase importer application")
    val cmd = createOptions().createCommandLine("knowledgebase-importer", args)

    val doidOwlLocation = cmd.getOptionValue(DOID_OWL_LOCATION)
    logger.info("Loading disease ontology from $doidOwlLocation")
    val diseaseOntology = DiseaseOntology(doidOwlLocation)

    logger.info("Reading knowledgebases...")
    val knowledgebases = readKnowledgebases(cmd, diseaseOntology)

    logger.info("Reading cancer type DOIDs from knowledgebases...")
    val cancerTypesDoids = knowledgebaseCancerDoids(knowledgebases, diseaseOntology)

    val outputDir = cmd.getOptionValue(OUTPUT_DIRECTORY)
    logger.info("Generating output files and writing to $outputDir")
    writeOutput(outputDir, knowledgebases, cancerTypesDoids)

    logger.info("Done")
}

private fun createOptions(): HmfOptions {
    val options = HmfOptions()
    options.add(RequiredInputFileOption(REFERENCE, "path to ref genome file"))
    options.add(RequiredCommandOption(TRANSVAR_LOCATION, "path to transvar location"))
    options.add(RequiredInputFileOption(DOID_OWL_LOCATION, "path to doid owl file"))
    options.add(RequiredInputFileOption(ONCO_ANNOTATED_LOCATION, "path to OncoKB annotated variants file"))
    options.add(RequiredInputFileOption(ONCO_ACTIONABLE_LOCATION, "path to OncoKB actionable variants file"))
    options.add(RequiredInputFileOption(CGI_VALIDATED_LOCATION, "path to cgi validated mutations file"))
    options.add(RequiredInputFileOption(CGI_BIOMARKERS_LOCATION, "path to cgi biomarkers file"))
    options.add(RequiredInputFileOption(CIVIC_VARIANTS_LOCATION, "path to civic variants file"))
    options.add(RequiredInputFileOption(CIVIC_EVIDENCE_LOCATION, "path to civic evidence file"))
    options.add(RequiredInputFileOption(COSMIC_FUSIONS_LOCATION, "path to cosmic fusions file"))
    options.add(RequiredInputFileOption(TREATMENT_TYPE_MAPPING_LOCATION, "path to treatment type mapping file"))
    options.add(RequiredInputOption(ICLUSION_ENDPOINT, "iclusion endpoint"))
    options.add(RequiredInputOption(ICLUSION_CLIENT_ID, "iclusion clientId"))
    options.add(RequiredInputOption(ICLUSION_CLIENT_SECRET, "iclusion client secret"))
    options.add(RequiredInputOption(ICLUSION_USER, "iclusion user"))
    options.add(RequiredInputOption(ICLUSION_PASSWORD, "iclusion password"))
    options.add(RequiredOutputOption(OUTPUT_DIRECTORY, "path to output directory"))
    return options
}

private fun readKnowledgebases(cmd: CommandLine, diseaseOntology: DiseaseOntology): List<Knowledgebase> {
    val reference = IndexedFastaSequenceFile(File(cmd.getOptionValue(REFERENCE)))
    val transvar = cmd.getOptionValue(TRANSVAR_LOCATION)
    val recordAnalyzer = RecordAnalyzer(transvar, reference)
    val treatmentTypeMap = CsvReader.readTSV<HmfDrug>(cmd.getOptionValue(TREATMENT_TYPE_MAPPING_LOCATION))
            .associateBy({ it.name.toLowerCase() }, { it.type })

    val oncoKb = OncoKb(cmd.getOptionValue(ONCO_ANNOTATED_LOCATION), cmd.getOptionValue(ONCO_ACTIONABLE_LOCATION), diseaseOntology,
            recordAnalyzer, treatmentTypeMap)
    val cgi = Cgi(cmd.getOptionValue(CGI_VALIDATED_LOCATION), cmd.getOptionValue(CGI_BIOMARKERS_LOCATION), diseaseOntology,
            recordAnalyzer, treatmentTypeMap)
    val civic = Civic(cmd.getOptionValue(CIVIC_VARIANTS_LOCATION), cmd.getOptionValue(CIVIC_EVIDENCE_LOCATION), diseaseOntology,
            recordAnalyzer, treatmentTypeMap)
    val cosmic = Cosmic(cmd.getOptionValue(COSMIC_FUSIONS_LOCATION))
    val iclusion = Iclusion(readIclusionStudies(cmd), diseaseOntology, recordAnalyzer)

//    return listOf(oncoKb, cgi, civic, cosmic, iclusion)
    return listOf(iclusion);
}

private fun readIclusionStudies(cmd: CommandLine): List<IclusionStudyDetails> {
    val iclusionEndpoint = cmd.getOptionValue(ICLUSION_ENDPOINT)
    logger.info("Connecting with iclusion API on $iclusionEndpoint")
    val iclusionApi = IclusionApiWrapper(iclusionEndpoint, cmd.getOptionValue(ICLUSION_CLIENT_ID),
            cmd.getOptionValue(ICLUSION_CLIENT_SECRET), cmd.getOptionValue(ICLUSION_USER), cmd.getOptionValue(ICLUSION_PASSWORD))
    logger.info("Reading iclusion study details...")
    val iclusionStudies = iclusionApi.studyDetails()
    iclusionApi.close()

    iclusionStudies.forEach {
        logger.info("iclusion study: $it")
    }
    logger.info("Queried and filtered ${iclusionStudies.size} studies from iclusion API")

    return iclusionStudies
}

private fun writeOutput(outputDir: String, knowledgebases: List<Knowledgebase>, cancerTypesDoids: List<CancerTypeDoidOutput>) {
    val dir = File(outputDir)
    if (!dir.exists()) dir.mkdirs()

    logger.info("Writing known variants to $outputDir")
    knowledgebases.filterNot { it.knownVariants.isEmpty() }.map { writeKnownVariants(it, outputDir) }

    val fusionPairLocation = "$outputDir${File.separator}knownFusionPairs.csv"
    logger.info("Writing known fusion genes to $fusionPairLocation")
    CsvWriter.writeCSV(knownFusionPairs(knowledgebases), fusionPairLocation)

    val promiscuousFiveGeneLocation = "$outputDir${File.separator}knownPromiscuousFive.csv"
    logger.info("Writing known promiscuous 5' genes to $promiscuousFiveGeneLocation")
    CsvWriter.writeCSV(knownPromiscuousFive(knowledgebases), promiscuousFiveGeneLocation)

    val promiscuousThreeGeneLocation = "$outputDir${File.separator}knownPromiscuousThree.csv"
    logger.info("Writing known promiscuous 3' genes to $promiscuousThreeGeneLocation")
    CsvWriter.writeCSV(knownPromiscuousThree(knowledgebases), promiscuousThreeGeneLocation)

    logger.info("Writing actionability files to $outputDir")
    CsvWriter.writeTSV(knowledgebases.flatMap { it.actionableVariants }, "$outputDir${File.separator}actionableVariants.tsv")
    CsvWriter.writeTSV(actionableFusionPairs(knowledgebases), "$outputDir${File.separator}actionableFusionPairs.tsv")
    CsvWriter.writeTSV(actionablePromiscuousFive(knowledgebases), "$outputDir${File.separator}actionablePromiscuousFive.tsv")
    CsvWriter.writeTSV(actionablePromiscuousThree(knowledgebases), "$outputDir${File.separator}actionablePromiscuousThree.tsv")
    CsvWriter.writeTSV(knowledgebases.flatMap { it.actionableCNVs }, "$outputDir${File.separator}actionableCNVs.tsv")
    CsvWriter.writeTSV(knowledgebases.flatMap { it.actionableRanges }, "$outputDir${File.separator}actionableRanges.tsv")

    val cancerTypeDOIDMappingLocation = "$outputDir${File.separator}knowledgebaseCancerTypes.tsv"
    logger.info("Writing cancer type DOID mapping to $cancerTypeDOIDMappingLocation")
    CsvWriter.writeTSV(cancerTypesDoids, cancerTypeDOIDMappingLocation)
}

private fun knowledgebaseCancerDoids(knowledgebases: List<Knowledgebase>, ontology: DiseaseOntology): List<CancerTypeDoidOutput> {
    val extraCancerTypeDoids = readExtraCancerTypeDoids().map {
        Pair(it.key, it.value.flatMap { doid -> ontology.findDoids(doid) }.toSet().sortedBy { it.value })
    }.toMap()
    val allCancerTypeDoids = knowledgebases.fold(mapOf<String, Set<Doid>>()) { map, it -> map + it.cancerTypes }
    return (allCancerTypeDoids + extraCancerTypeDoids).entries.map { CancerTypeDoidOutput(it.key, it.value.joinToString(";")) }
}

private fun readExtraCancerTypeDoids(): Map<String, Set<Doid>> {
    return readCSVRecords(object {}.javaClass.getResourceAsStream("/knowledgebase_disease_doids.csv")) {
        Pair(it["cancerType"], it["doids"].orEmpty().split(";").filterNot { it.isBlank() }.map { Doid(it.trim()) }.toSet())
    }.toMap()
}

private fun writeKnownVariants(knowledgebase: Knowledgebase, outputDirectory: String) {
    CsvWriter.writeTSV(knowledgebase.knownVariants.distinct(), "$outputDirectory${File.separator}${knowledgebase.source}KnownVariants.tsv")
}
