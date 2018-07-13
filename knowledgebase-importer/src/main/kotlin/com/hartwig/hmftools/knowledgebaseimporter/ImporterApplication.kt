package com.hartwig.hmftools.knowledgebaseimporter

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
import com.hartwig.hmftools.knowledgebaseimporter.dao.EnsemblGeneDAO
import com.hartwig.hmftools.knowledgebaseimporter.diseaseOntology.DiseaseOntology
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.RecordAnalyzer
import com.hartwig.hmftools.knowledgebaseimporter.oncoKb.OncoKb
import com.hartwig.hmftools.knowledgebaseimporter.output.CancerTypeDoidOutput
import com.hartwig.hmftools.knowledgebaseimporter.output.HmfDrug
import htsjdk.samtools.reference.IndexedFastaSequenceFile
import org.apache.commons.cli.CommandLine
import java.io.File

fun main(args: Array<String>) {
    val cmd = createOptions().createCommandLine("knowledgebase-importer", args)
    val outputDir = cmd.getOptionValue(OUTPUT_DIRECTORY)
    val diseaseOntology = DiseaseOntology(cmd.getOptionValue(DOID_OWL_LOCATION))
    val knowledgebases = readKnowledgebases(cmd, diseaseOntology)
    val cancerTypesDoids = knowledgebaseCancerDoids(knowledgebases, diseaseOntology)
    writeOutput(outputDir, knowledgebases, cancerTypesDoids)
}

private fun createOptions(): HmfOptions {
    val options = HmfOptions()
    options.add(RequiredInputFileOption(REFERENCE, "path to ref genome file"))
    options.add(RequiredCommandOption(TRANSVAR_LOCATION, "path to transvar location"))
    options.add(RequiredInputFileOption(DOID_OWL_LOCATION, "path to doid owl file"))
    options.add(RequiredInputFileOption(ONCO_ANNOTATED_LOCATION, "path to oncoKb annotated variants file"))
    options.add(RequiredInputFileOption(ONCO_ACTIONABLE_LOCATION, "path to oncoKB actionable variants file"))
    options.add(RequiredInputFileOption(CGI_VALIDATED_LOCATION, "path to cgi validated mutations file"))
    options.add(RequiredInputFileOption(CGI_BIOMARKERS_LOCATION, "path to cgi biomarkers file"))
    options.add(RequiredInputFileOption(CIVIC_VARIANTS_LOCATION, "path to civic variants file"))
    options.add(RequiredInputFileOption(CIVIC_EVIDENCE_LOCATION, "path to civic evidence file"))
    options.add(RequiredInputFileOption(COSMIC_FUSIONS_LOCATION, "path to cosmic fusions file"))
    options.add(RequiredInputFileOption(TREATMENT_TYPE_MAPPING_LOCATION, "path to treatment type mapping file"))
    options.add(RequiredOutputOption(OUTPUT_DIRECTORY, "path to output directory"))
    options.add(RequiredInputOption(ENSEMBL_DB, "ensembl db url"))
    options.add(RequiredInputOption(DB_USER, "db user"))
    options.add(RequiredInputOption(DB_PASSWORD, "db password"))
    return options
}

private fun readKnowledgebases(cmd: CommandLine, diseaseOntology: DiseaseOntology): List<Knowledgebase> {
    val reference = IndexedFastaSequenceFile(File(cmd.getOptionValue(REFERENCE)))
    val ensemblJdbcUrl = "jdbc:${cmd.getOptionValue(ENSEMBL_DB)}"
    val ensemblGeneDAO = EnsemblGeneDAO(ensemblJdbcUrl, cmd.getOptionValue(DB_USER), cmd.getOptionValue(DB_PASSWORD))
    val transvar = cmd.getOptionValue(TRANSVAR_LOCATION)
    val recordAnalyzer = RecordAnalyzer(transvar, reference, ensemblGeneDAO)
    val treatmentTypeMap = CsvReader.readTSV<HmfDrug>(cmd.getOptionValue(TREATMENT_TYPE_MAPPING_LOCATION))
            .associateBy({ it.name.toLowerCase() }, { it.type })

    val oncoKb = OncoKb(cmd.getOptionValue(ONCO_ANNOTATED_LOCATION), cmd.getOptionValue(ONCO_ACTIONABLE_LOCATION), diseaseOntology,
                        recordAnalyzer, treatmentTypeMap)
    val cgi = Cgi(cmd.getOptionValue(CGI_VALIDATED_LOCATION), cmd.getOptionValue(CGI_BIOMARKERS_LOCATION), diseaseOntology,
                  recordAnalyzer, treatmentTypeMap)
    val civic = Civic(cmd.getOptionValue(CIVIC_VARIANTS_LOCATION), cmd.getOptionValue(CIVIC_EVIDENCE_LOCATION), diseaseOntology,
                      recordAnalyzer, treatmentTypeMap)
    val cosmic = Cosmic(cmd.getOptionValue(COSMIC_FUSIONS_LOCATION))
    return listOf(oncoKb, cgi, civic, cosmic)
}

private fun writeOutput(outputDir: String, knowledgebases: List<Knowledgebase>, cancerTypesDoids: List<CancerTypeDoidOutput>) {
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
    CsvWriter.writeTSV(cancerTypesDoids, "$outputDir${File.separator}knowledgebaseCancerTypes.tsv")
}

private fun knowledgebaseCancerDoids(knowledgebases: List<Knowledgebase>, ontology: DiseaseOntology): List<CancerTypeDoidOutput> {
    val extraCancerTypeDoids = readExtraCancerTypeDoids().map {
        Pair(it.key, it.value.flatMap { doid -> ontology.findDoidsForDoid(doid) }.toSet().sorted())
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
