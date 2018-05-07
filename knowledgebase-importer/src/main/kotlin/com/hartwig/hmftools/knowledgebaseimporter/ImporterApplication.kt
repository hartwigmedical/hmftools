package com.hartwig.hmftools.knowledgebaseimporter

import com.hartwig.hmftools.extensions.cli.createCommandLine
import com.hartwig.hmftools.knowledgebaseimporter.cgi.Cgi
import com.hartwig.hmftools.knowledgebaseimporter.civic.Civic
import com.hartwig.hmftools.knowledgebaseimporter.cosmic.Cosmic
import com.hartwig.hmftools.knowledgebaseimporter.diseaseOntology.DiseaseOntology
import com.hartwig.hmftools.knowledgebaseimporter.oncoKb.OncoKb
import com.hartwig.hmftools.knowledgebaseimporter.output.KnownVariantOutput
import htsjdk.samtools.reference.IndexedFastaSequenceFile
import org.apache.commons.cli.Option
import org.apache.commons.cli.Options
import org.apache.commons.csv.CSVFormat
import org.apache.commons.csv.CSVPrinter
import java.io.File
import java.io.FileWriter

fun main(args: Array<String>) {
    val cmd = createOptions().createCommandLine("knowledgebase-importer", args)
    val reference = IndexedFastaSequenceFile(File(cmd.getOptionValue(REFERENCE)))
    val diseaseOntology = DiseaseOntology(cmd.getOptionValue(DOID_OWL_LOCATION))
    val transvar = cmd.getOptionValue(TRANSVAR_LOCATION)
    val oncoKb = OncoKb(cmd.getOptionValue(ONCO_ANNOTATED_LOCATION), cmd.getOptionValue(ONCO_ACTIONABLE_LOCATION), transvar,
                        diseaseOntology, reference)
    val cgi = Cgi(cmd.getOptionValue(CGI_VALIDATED_LOCATION), cmd.getOptionValue(CGI_BIOMARKERS_LOCATION), transvar, diseaseOntology,
                  reference)
    val civic = Civic(cmd.getOptionValue(CIVIC_VARIANTS_LOCATION), cmd.getOptionValue(CIVIC_EVIDENCE_LOCATION), transvar, diseaseOntology,
                      reference)
    val cosmic = Cosmic("cosmic_gene_fusions.csv")
    val knowledgebases = listOf(oncoKb, cgi, civic, cosmic)
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

private fun writeKnownVariants(knowledgebase: Knowledgebase, outputDirectory: String) {
    val format = CSVFormat.TDF.withHeader(*KnownVariantOutput.header.toTypedArray()).withNullString("")
    val printer = CSVPrinter(FileWriter("$outputDirectory${File.separator}${knowledgebase.source}_known_variants"), format)
    printer.printRecords(knowledgebase.knownVariants.distinct().map { it.record })
    printer.close()
}
