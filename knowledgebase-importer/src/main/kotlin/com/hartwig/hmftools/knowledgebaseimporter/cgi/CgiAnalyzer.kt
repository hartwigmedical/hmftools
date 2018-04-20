package com.hartwig.hmftools.knowledgebaseimporter.cgi

import com.hartwig.hmftools.knowledgebaseimporter.output.KnownVariantOutput
import com.hartwig.hmftools.knowledgebaseimporter.transvar.TransvarProteinAnalyzer
import com.hartwig.hmftools.knowledgebaseimporter.transvar.extract
import com.hartwig.hmftools.knowledgebaseimporter.transvar.extractChromosome
import com.hartwig.hmftools.knowledgebaseimporter.transvar.extractVariants
import htsjdk.samtools.reference.IndexedFastaSequenceFile
import org.apache.commons.csv.CSVFormat
import org.apache.commons.csv.CSVParser
import java.io.File
import java.nio.charset.Charset

fun analyzeCgi(transvarLocation: String, fileLocation: String, reference: IndexedFastaSequenceFile): List<KnownVariantOutput> {
    val parser = CSVParser.parse(File(fileLocation), Charset.defaultCharset(), CSVFormat.TDF.withFirstRecordAsHeader().withNullString(""))
    val analyzer = TransvarProteinAnalyzer(transvarLocation)
    val records = parser.asSequence().map { CgiRecord(it) }.toList().filter { it.context == "somatic" }
    val transvarOutput = analyzer.analyze(records.map { Pair(it.transcript, it.impact) })
    return records.zip(transvarOutput)
            .flatMap { (cgiRecord, transvarOutput) ->
                val cgiVariantGdnas = cgiRecord.gdna.split("__").filterNot { it.isBlank() }
                val inferredVariants = extractVariants(transvarOutput, reference)
                val chromosome = extractChromosome(cgiVariantGdnas.first())
                val cgiVariants = extract(chromosome, cgiVariantGdnas, reference)
                val cgiVariantOutputs = cgiVariants.map { KnownVariantOutput(cgiRecord.gene, cgiRecord.transcript, "CGI", it) }
                val inferredVariantOutputs = inferredVariants.filterNot { variant -> cgiVariants.any { it == variant } }.map {
                    KnownVariantOutput(cgiRecord.gene, cgiRecord.transcript, "", it)
                }
                cgiVariantOutputs + inferredVariantOutputs
            }
}
