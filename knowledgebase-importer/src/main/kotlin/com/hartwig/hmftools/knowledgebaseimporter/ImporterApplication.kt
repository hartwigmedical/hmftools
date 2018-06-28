package com.hartwig.hmftools.knowledgebaseimporter

import com.hartwig.hmftools.extensions.cli.createCommandLine
import com.hartwig.hmftools.extensions.csv.CsvReader
import com.hartwig.hmftools.extensions.csv.CsvWriter
import com.hartwig.hmftools.knowledgebaseimporter.cgi.Cgi
import com.hartwig.hmftools.knowledgebaseimporter.civic.Civic
import com.hartwig.hmftools.knowledgebaseimporter.cosmic.Cosmic
import com.hartwig.hmftools.knowledgebaseimporter.dao.EnsemblGeneDAO
import com.hartwig.hmftools.knowledgebaseimporter.diseaseOntology.DiseaseOntology
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.RecordAnalyzer
import com.hartwig.hmftools.knowledgebaseimporter.oncoKb.OncoKb
import com.hartwig.hmftools.knowledgebaseimporter.output.CancerTypeDoidOutput
import com.hartwig.hmftools.knowledgebaseimporter.output.HmfDrug
import htsjdk.samtools.reference.IndexedFastaSequenceFile
import org.apache.commons.cli.Option
import org.apache.commons.cli.Options
import java.io.File

fun main(args: Array<String>) {
    val cmd = createOptions().createCommandLine("knowledgebase-importer", args)
    val outputDir = cmd.getOptionValue(OUTPUT_DIRECTORY)
    val reference = IndexedFastaSequenceFile(File(cmd.getOptionValue(REFERENCE)))
    val diseaseOntology = DiseaseOntology(cmd.getOptionValue(DOID_OWL_LOCATION))
    val treatmentTypeMap = CsvReader.readTSV<HmfDrug>("treatmentTypeMapping.tsv")
            .associateBy({ it.name.toLowerCase() }, { it.type })
    val ensemblGeneDAO = EnsemblGeneDAO(cmd.getOptionValue(ENSEMBL_URL), cmd.getOptionValue(DB_USER), cmd.getOptionValue(DB_PASSWORD))
    val transvar = cmd.getOptionValue(TRANSVAR_LOCATION)
    val recordAnalyzer = RecordAnalyzer(transvar, reference, ensemblGeneDAO)

    val oncoKb = OncoKb(cmd.getOptionValue(ONCO_ANNOTATED_LOCATION), cmd.getOptionValue(ONCO_ACTIONABLE_LOCATION), diseaseOntology,
                        recordAnalyzer, treatmentTypeMap)
    val cgi = Cgi(cmd.getOptionValue(CGI_VALIDATED_LOCATION), cmd.getOptionValue(CGI_BIOMARKERS_LOCATION), diseaseOntology,
                  recordAnalyzer, treatmentTypeMap)
    val civic = Civic(cmd.getOptionValue(CIVIC_VARIANTS_LOCATION), cmd.getOptionValue(CIVIC_EVIDENCE_LOCATION), diseaseOntology,
                      recordAnalyzer, treatmentTypeMap)
    val cosmic = Cosmic("cosmic_gene_fusions.csv")
    val knowledgebases = listOf(oncoKb, cgi, civic, cosmic)
    val cancerTypesDoids = knowledgebaseCancerDoids(knowledgebases, diseaseOntology)

    knowledgebases.filterNot { it.knownVariants.isEmpty() }.map { writeKnownVariants(it, outputDir) }
    CsvWriter.writeCSV(knownFusionPairs(knowledgebases), "$outputDir${File.separator}knownFusionPairs.csv")
    CsvWriter.writeCSV(knownPromiscuousFive(knowledgebases), "$outputDir${File.separator}knownPromiscuousFive.csv")
    CsvWriter.writeCSV(knownPromiscuousThree(knowledgebases), "$outputDir${File.separator}knownPromiscuousFive.csv")

    CsvWriter.writeTSV(knowledgebases.flatMap { it.actionableVariants }, "$outputDir${File.separator}actionableVariants.tsv")
    CsvWriter.writeTSV(actionableFusionPairs(knowledgebases), "$outputDir${File.separator}actionableFusionPairs.tsv")
    CsvWriter.writeTSV(actionablePromiscuousFive(knowledgebases), "$outputDir${File.separator}actionablePromiscuousFive.tsv")
    CsvWriter.writeTSV(actionablePromiscuousThree(knowledgebases), "$outputDir${File.separator}actionablePromiscuousThree.tsv")
    CsvWriter.writeTSV(knowledgebases.flatMap { it.actionableCNVs }, "$outputDir${File.separator}actionableCNVs.tsv")
    CsvWriter.writeTSV(knowledgebases.flatMap { it.actionableRanges }, "$outputDir${File.separator}actionableRanges.tsv")
    CsvWriter.writeTSV(knowledgebases.flatMap { it.actionableGenes }, "$outputDir${File.separator}actionableGenes.tsv")
    CsvWriter.writeTSV(cancerTypesDoids, "$outputDir${File.separator}knowledgebaseCancerTypes.tsv")
}

private fun createOptions(): Options {
    val options = Options()
    options.addOption(Option.builder(REFERENCE).required().hasArg().desc("path to ref genome file").build())
    options.addOption(Option.builder(TRANSVAR_LOCATION).required().hasArg().desc("path to transvar location").build())
    options.addOption(Option.builder(DOID_OWL_LOCATION).required().hasArg().desc("path to doid owl file").build())
    options.addOption(Option.builder(ONCO_ANNOTATED_LOCATION).required().hasArg().desc("path to oncoKb annotated variants file").build())
    options.addOption(Option.builder(ONCO_ACTIONABLE_LOCATION).required().hasArg().desc("path to oncoKB actionable variants file").build())
    options.addOption(Option.builder(CGI_VALIDATED_LOCATION).required().hasArg().desc("path to cgi validated mutations file").build())
    options.addOption(Option.builder(CGI_BIOMARKERS_LOCATION).required().hasArg().desc("path to cgi biomarkers file").build())
    options.addOption(Option.builder(CIVIC_VARIANTS_LOCATION).required().hasArg().desc("path to civic variants file").build())
    options.addOption(Option.builder(CIVIC_EVIDENCE_LOCATION).required().hasArg().desc("path to civic evidence file").build())
    options.addOption(Option.builder(COSMIC_FUSIONS_LOCATION).required().hasArg().desc("path to cosmic fusions file").build())
    options.addOption(Option.builder(OUTPUT_DIRECTORY).required().hasArg().desc("path to output directory").build())
    options.addOption(Option.builder(ENSEMBL_URL).required().hasArg().desc("ensembl db url").build())
    options.addOption(Option.builder(DB_USER).required().hasArg().desc("db user").build())
    options.addOption(Option.builder(DB_PASSWORD).required().hasArg().desc("db password").build())
    return options
}

private fun knowledgebaseCancerDoids(knowledgebases: List<Knowledgebase>, ontology: DiseaseOntology): List<CancerTypeDoidOutput> {
    val extraCancerTypeDoids = readExtraCancerTypeDoids().map {
        Pair(it.key, it.value.flatMap { doid -> ontology.findDoidsForDoid(doid) }.toSet())
    }.toMap()
    val allCancerTypeDoids = knowledgebases.fold(mapOf<String, Set<String>>(), { map, it -> map + it.cancerTypes })
    return (allCancerTypeDoids + extraCancerTypeDoids).entries.map { CancerTypeDoidOutput(it.key, it.value.joinToString(";")) }
}

private fun readExtraCancerTypeDoids(): Map<String, Set<String>> {
    return readCSVRecords(object {}.javaClass.getResourceAsStream("/knowledgebase_disease_doids.csv")) {
        Pair(it["cancerType"], it["doids"].orEmpty().split(";").filterNot { it.isBlank() }.map { it.trim() }.toSet())
    }.toMap()
}

private fun writeKnownVariants(knowledgebase: Knowledgebase, outputDirectory: String) {
    CsvWriter.writeTSV(knowledgebase.knownVariants.distinct(), "$outputDirectory${File.separator}${knowledgebase.source}KnownVariants.tsv")
}
