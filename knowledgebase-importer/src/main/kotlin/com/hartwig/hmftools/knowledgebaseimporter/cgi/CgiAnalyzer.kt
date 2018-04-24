package com.hartwig.hmftools.knowledgebaseimporter.cgi

import com.hartwig.hmftools.common.variant.SomaticVariant
import com.hartwig.hmftools.knowledgebaseimporter.output.Actionability
import com.hartwig.hmftools.knowledgebaseimporter.output.ActionableVariantOutput
import com.hartwig.hmftools.knowledgebaseimporter.output.KnownVariantOutput
import com.hartwig.hmftools.knowledgebaseimporter.transvar.*
import com.hartwig.hmftools.knowledgebaseimporter.transvar.annotations.CDnaAnnotation
import com.hartwig.hmftools.knowledgebaseimporter.transvar.annotations.ProteinAnnotation
import htsjdk.samtools.reference.IndexedFastaSequenceFile
import org.apache.commons.csv.CSVFormat
import org.apache.commons.csv.CSVParser
import org.apache.commons.csv.CSVRecord
import java.io.File
import java.nio.charset.Charset

fun analyzeCgi(transvarLocation: String, fileLocation: String, reference: IndexedFastaSequenceFile): List<KnownVariantOutput> {
    val somaticRecords = readCgiRecords(fileLocation) { CgiKnownVariantRecord(it) }.filter { it.context == "somatic" }
    val transvarOutput = annotateCgiRecords(transvarLocation, somaticRecords)
    return somaticRecords.zip(transvarOutput)
            .flatMap { (cgiRecord, transvarOutput) ->
                val cgiVariants = extractCgiVariants(cgiRecord.gdna, reference)
                val inferredVariants = extractVariants(transvarOutput, reference).filterNot { variant -> cgiVariants.any { it == variant } }
                val cgiVariantOutputs = cgiVariants.map { KnownVariantOutput(cgiRecord.gene, cgiRecord.transcript, "CGI", it) }
                val inferredVariantOutputs = inferredVariants.map { KnownVariantOutput(cgiRecord.gene, cgiRecord.transcript, "", it) }
                cgiVariantOutputs + inferredVariantOutputs
            }
}

fun <T> readCgiRecords(fileLocation: String, readRecord: (CSVRecord) -> T): List<T> {
    val parser = CSVParser.parse(File(fileLocation), Charset.defaultCharset(), CSVFormat.TDF.withFirstRecordAsHeader().withNullString(""))
    return parser.asSequence().map { readRecord(it) }.toList()
}

fun annotateCgiRecords(transvarLocation: String, records: List<CgiKnownVariantRecord>): List<TransvarOutput> {
    val analyzer = TransvarProteinAnalyzer(transvarLocation)
    return analyzer.analyze(records.map { ProteinAnnotation(it.transcript, it.impact) })
}

fun annotateCgiBiomarkers(transvarLocation: String, records: List<CgiBiomarkersRecord>): List<Pair<TransvarOutput, TransvarOutput>> {
    val proteinAnalyzer = TransvarProteinAnalyzer(transvarLocation)
    val cdnaAnalyzer = TransvarCdnaAnalyzer(transvarLocation)
    val proteinOutput = proteinAnalyzer.analyze(records.map { ProteinAnnotation(it.transcript, it.protein) })
    val cdnaOutput = cdnaAnalyzer.analyze(records.map { CDnaAnnotation(it.transcript, it.cdna) })
    return proteinOutput.zip(cdnaOutput)
}

fun extractCgiVariants(gdna: String, reference: IndexedFastaSequenceFile): List<SomaticVariant> {
    val cgiVariantGdnas = gdna.split("__").filterNot { it.isBlank() }
    val chromosome = extractChromosome(cgiVariantGdnas.first())
    return extract(chromosome, cgiVariantGdnas, reference)
}

fun analyzeCgiActionable(transvarLocation: String, fileLocation: String,
                         reference: IndexedFastaSequenceFile): List<ActionableVariantOutput> {
    val records = readCgiRecords(fileLocation) { CgiBiomarkersRecord(it) }
    val transvarOutput = annotateCgiBiomarkers(transvarLocation, records)
    return records.zip(transvarOutput)
            .flatMap { (cgiRecord, transvarOutput) ->
                val cgiVariants = extractCgiVariants(cgiRecord.gdna, reference)
                val inferredProteinVariants = extractVariants(transvarOutput.first, reference)
                        .filterNot { variant -> cgiVariants.any { it == variant } }
                val inferredCdnaVariants = extractVariants(transvarOutput.second, reference)
                        .filterNot { variant -> cgiVariants.any { it == variant } || inferredProteinVariants.any { it == variant } }
                actionableVariantOutput(cgiRecord, cgiVariants + inferredProteinVariants + inferredCdnaVariants)
            }
}

fun actionableVariantOutput(cgiRecord: CgiBiomarkersRecord, somaticVariants: List<SomaticVariant>): List<ActionableVariantOutput> {
    return cgiRecord.cancerTypes.flatMap { cancerType ->
        somaticVariants.map { somaticVariant ->
            ActionableVariantOutput(cgiRecord.gene, somaticVariant,
                                    Actionability("cgi", cancerType, cgiRecord.drug, cgiRecord.level, cgiRecord.association, ""))
        }
    }
}
