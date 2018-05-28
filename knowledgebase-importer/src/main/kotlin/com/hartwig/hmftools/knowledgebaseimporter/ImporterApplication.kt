package com.hartwig.hmftools.knowledgebaseimporter

import com.hartwig.hmftools.extensions.cli.createCommandLine
import com.hartwig.hmftools.knowledgebaseimporter.cgi.Cgi
import com.hartwig.hmftools.knowledgebaseimporter.civic.Civic
import com.hartwig.hmftools.knowledgebaseimporter.cosmic.Cosmic
import com.hartwig.hmftools.knowledgebaseimporter.diseaseOntology.DiseaseOntology
import com.hartwig.hmftools.knowledgebaseimporter.oncoKb.OncoKb
import com.hartwig.hmftools.knowledgebaseimporter.output.*
import htsjdk.samtools.reference.IndexedFastaSequenceFile
import org.apache.commons.cli.Option
import org.apache.commons.cli.Options
import org.apache.commons.csv.CSVFormat
import org.apache.commons.csv.CSVPrinter
import java.io.File
import java.io.FileWriter

fun main(args: Array<String>) {
    val cmd = createOptions().createCommandLine("knowledgebase-importer", args)
    val outputDir = cmd.getOptionValue(OUTPUT_DIRECTORY)
    val reference = IndexedFastaSequenceFile(File(cmd.getOptionValue(REFERENCE)))
    val diseaseOntology = DiseaseOntology(cmd.getOptionValue(DOID_OWL_LOCATION))
    val treatmentTypeMap = readTreatmentTypes("treatmentTypeMapping.tsv")

    val transvar = cmd.getOptionValue(TRANSVAR_LOCATION)
    val oncoKb = OncoKb(cmd.getOptionValue(ONCO_ANNOTATED_LOCATION), cmd.getOptionValue(ONCO_ACTIONABLE_LOCATION), transvar,
                        diseaseOntology, reference, treatmentTypeMap)
    val cgi = Cgi(cmd.getOptionValue(CGI_VALIDATED_LOCATION), cmd.getOptionValue(CGI_BIOMARKERS_LOCATION), transvar, diseaseOntology,
                  reference, treatmentTypeMap)
    val civic = Civic(cmd.getOptionValue(CIVIC_VARIANTS_LOCATION), cmd.getOptionValue(CIVIC_EVIDENCE_LOCATION), transvar, diseaseOntology,
                      reference, treatmentTypeMap)
    val cosmic = Cosmic("cosmic_gene_fusions.csv")
    val knowledgebases = listOf(oncoKb, cgi, civic, cosmic)
    val cancerTypesDoids = knowledgebaseCancerDoids(knowledgebases, diseaseOntology)

    knowledgebases.filterNot { it.knownVariants.isEmpty() }.map { writeKnownVariants(it, outputDir) }
    writeKnownFusionPairs(knownFusionPairs(knowledgebases), "$outputDir${File.separator}knownFusionPairs.csv")
    writePromiscuousGenes(knownPromiscuousFive(knowledgebases), "$outputDir${File.separator}knownPromiscuousFive.csv")
    writePromiscuousGenes(knownPromiscuousThree(knowledgebases), "$outputDir${File.separator}knownPromiscuousFive.csv")

    writeActionableVariants(knowledgebases.flatMap { it.actionableVariants }, "$outputDir${File.separator}actionableVariants.tsv")
    writeActionableFusionPairs(actionableFusionPairs(knowledgebases), "$outputDir${File.separator}actionableFusionPairs.tsv")
    writeActionablePromiscuousGenes(actionablePromiscuousFive(knowledgebases), "$outputDir${File.separator}actionablePromiscuousFive.tsv")
    writeActionablePromiscuousGenes(actionablePromiscuousThree(knowledgebases), "$outputDir${File.separator}actionablePromiscuousThree.tsv")
    writeActionableCnvs(knowledgebases.flatMap { it.actionableCNVs }, "$outputDir${File.separator}actionableCNVs.tsv")
    writeCancerDoids(cancerTypesDoids, "$outputDir${File.separator}knowledgebaseCancerTypes.tsv")
    //    writeTreatmentTypes(bootstrapTreatmentTypeMapping(cgi), "$outputDir${File.separator}treatmentTypeMapping.tsv")
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

private fun bootstrapTreatmentTypeMapping(cgi: Cgi): List<HmfDrug> {
    return cgi.actionableKbRecords.flatMap { it.cgiDrugs }.groupBy { it.drugName }
            .mapValues { (key, value) ->
                val drugFamilies = value.flatMap { it.drugTypes }.toSet()
                if (drugFamilies.size > 1) drugFamilies.filterNot { it == key } else drugFamilies
            }.map { HmfDrug(it.key, it.value.joinToString(";")) }
}

private fun writeKnownVariants(knowledgebase: Knowledgebase, outputDirectory: String) {
    val format = CSVFormat.TDF.withHeader(*KnownVariantOutput.header.toTypedArray()).withNullString("")
    val printer = CSVPrinter(FileWriter("$outputDirectory${File.separator}${knowledgebase.source}KnownVariants.tsv"), format)
    printer.printRecords(knowledgebase.knownVariants.distinct().map { it.record })
    printer.close()
}

private fun writeKnownFusionPairs(fusions: List<FusionPair>, location: String) {
    val format = CSVFormat.DEFAULT.withHeader(*FusionPair.header.toTypedArray()).withNullString("")
    val printer = CSVPrinter(FileWriter(location), format)
    printer.printRecords(fusions.distinct().sortedBy { it.fiveGene + it.threeGene }.map { it.record })
    printer.close()
}

private fun writePromiscuousGenes(fusions: List<PromiscuousGene>, location: String) {
    val format = CSVFormat.DEFAULT.withHeader(*PromiscuousGene.header.toTypedArray()).withNullString("")
    val printer = CSVPrinter(FileWriter(location), format)
    printer.printRecords(fusions.distinct().sortedBy { it.gene }.map { it.record })
    printer.close()
}

private fun writeActionableVariants(outputVariants: List<ActionableVariantOutput>, location: String) {
    val format = CSVFormat.TDF.withHeader(*ActionableVariantOutput.header.toTypedArray()).withNullString("")
    val printer = CSVPrinter(FileWriter(location), format)
    printer.printRecords(outputVariants.distinct().map { it.record })
    printer.close()
}

private fun writeActionableCnvs(cnvs: List<ActionableCNVOutput>, location: String) {
    val format = CSVFormat.TDF.withHeader(*ActionableCNVOutput.header.toTypedArray()).withNullString("")
    val printer = CSVPrinter(FileWriter(location), format)
    printer.printRecords(cnvs.distinct().map { it.record })
    printer.close()
}

private fun writeActionableFusionPairs(fusions: List<ActionableFusionOutput>, location: String) {
    val format = CSVFormat.TDF.withHeader(*ActionableFusionOutput.fusionPairHeader.toTypedArray()).withNullString("")
    val printer = CSVPrinter(FileWriter(location), format)
    printer.printRecords(fusions.map { it.record })
    printer.close()
}

private fun writeActionablePromiscuousGenes(fusions: List<ActionableFusionOutput>, location: String) {
    val format = CSVFormat.TDF.withHeader(*ActionableFusionOutput.promiscuousGeneHeader.toTypedArray()).withNullString("")
    val printer = CSVPrinter(FileWriter(location), format)
    printer.printRecords(fusions.map { it.record })
    printer.close()
}

private fun writeCancerDoids(cancerTypes: List<CancerTypeDoidOutput>, location: String) {
    val format = CSVFormat.TDF.withHeader(*CancerTypeDoidOutput.header.toTypedArray()).withNullString("")
    val printer = CSVPrinter(FileWriter(location), format)
    printer.printRecords(cancerTypes.map { it.record })
    printer.close()
}

private fun writeTreatmentTypes(output: List<HmfDrug>, location: String) {
    val format = CSVFormat.TDF.withHeader(*HmfDrug.header.toTypedArray()).withNullString("")
    val printer = CSVPrinter(FileWriter(location), format)
    printer.printRecords(output.distinct().map { it.record })
    printer.close()
}

private fun readTreatmentTypes(fileLocation: String): Map<String, String> {
    return readCSVRecords(fileLocation) {
        HmfDrug(it["name"].toLowerCase(), it["type"])
    }.associateBy({ it.name }, { it.type })
}
