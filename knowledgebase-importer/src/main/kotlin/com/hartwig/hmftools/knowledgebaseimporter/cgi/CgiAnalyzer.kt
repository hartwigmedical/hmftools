package com.hartwig.hmftools.knowledgebaseimporter.cgi

import com.hartwig.hmftools.common.variant.SomaticVariant
import com.hartwig.hmftools.knowledgebaseimporter.output.KnownVariantOutput
import com.hartwig.hmftools.knowledgebaseimporter.transvar.*
import com.hartwig.hmftools.knowledgebaseimporter.transvar.annotations.ProteinAnnotation
import htsjdk.samtools.reference.IndexedFastaSequenceFile
import org.apache.commons.csv.CSVFormat
import org.apache.commons.csv.CSVParser
import java.io.File
import java.nio.charset.Charset

fun analyzeCgi(transvarLocation: String, fileLocation: String, reference: IndexedFastaSequenceFile): List<KnownVariantOutput> {
    val records = readCgiRecords(fileLocation)
    val transvarOutput = annotateCgiRecords(transvarLocation, records)
    return records.zip(transvarOutput)
            .flatMap { (cgiRecord, transvarOutput) ->
                val cgiVariants = extractCgiVariants(cgiRecord, reference)
                val inferredVariants = extractVariants(transvarOutput, reference).filterNot { variant -> cgiVariants.any { it == variant } }
                val cgiVariantOutputs = cgiVariants.map { KnownVariantOutput(cgiRecord.gene, cgiRecord.transcript, "CGI", it) }
                val inferredVariantOutputs = inferredVariants.map { KnownVariantOutput(cgiRecord.gene, cgiRecord.transcript, "", it) }
                cgiVariantOutputs + inferredVariantOutputs
            }
}

fun readCgiRecords(fileLocation: String): List<CgiRecord> {
    val parser = CSVParser.parse(File(fileLocation), Charset.defaultCharset(), CSVFormat.TDF.withFirstRecordAsHeader().withNullString(""))
    return parser.asSequence().map { CgiRecord(it) }.toList().filter { it.context == "somatic" }
}

fun annotateCgiRecords(transvarLocation: String, records: List<CgiRecord>): List<TransvarOutput> {
    val analyzer = TransvarProteinAnalyzer(transvarLocation)
    return analyzer.analyze(records.map { ProteinAnnotation(it.transcript, it.impact) })
}

fun extractCgiVariants(cgiRecord: CgiRecord, reference: IndexedFastaSequenceFile): List<SomaticVariant> {
    val cgiVariantGdnas = cgiRecord.gdna.split("__").filterNot { it.isBlank() }
    val chromosome = extractChromosome(cgiVariantGdnas.first())
    return extract(chromosome, cgiVariantGdnas, reference)
}
