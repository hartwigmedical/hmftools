package com.hartwig.hmftools.knowledgebaseimporter.oncoKb

import com.hartwig.hmftools.knowledgebaseimporter.output.Actionability
import com.hartwig.hmftools.knowledgebaseimporter.output.ActionableVariantOutput
import com.hartwig.hmftools.knowledgebaseimporter.output.KnownVariantOutput
import com.hartwig.hmftools.knowledgebaseimporter.transvar.TransvarProteinAnalyzer
import com.hartwig.hmftools.knowledgebaseimporter.transvar.annotations.ProteinAnnotation
import com.hartwig.hmftools.knowledgebaseimporter.transvar.extractVariants
import htsjdk.samtools.reference.IndexedFastaSequenceFile
import org.apache.commons.csv.CSVFormat
import org.apache.commons.csv.CSVParser
import org.apache.commons.csv.CSVRecord
import java.io.File
import java.nio.charset.Charset

fun analyzeOncoKb(transvarLocation: String, fileLocation: String, reference: IndexedFastaSequenceFile): List<KnownVariantOutput> {
    val records = readOncoKbRecords(fileLocation) { OncoAnnotatedVariantRecord(it) }
    val analyzer = TransvarProteinAnalyzer(transvarLocation)
    val transvarOutput = analyzer.analyze(records.map { ProteinAnnotation(it.transcript, it.alteration) })
    return records.zip(transvarOutput)
            .flatMap { (oncoRecord, transvarOutput) ->
                extractVariants(transvarOutput, reference)
                        .map { KnownVariantOutput(oncoRecord.gene, oncoRecord.transcript, oncoRecord.oncogenicity, it) }
            }
}

fun analyzeOncoKbActionable(transvarLocation: String, fileLocation: String,
                            reference: IndexedFastaSequenceFile): List<ActionableVariantOutput> {
    val records = readOncoKbRecords(fileLocation) { OncoActionableVariantRecord(it) }
    return analyzePointMutations(transvarLocation, records, reference)
}

fun <T> readOncoKbRecords(fileLocation: String, recordParser: (CSVRecord) -> T): List<T> {
    val parser = CSVParser.parse(File(fileLocation), Charset.defaultCharset(), CSVFormat.TDF.withFirstRecordAsHeader().withNullString(""))
    return parser.asSequence().map { recordParser(it) }.toList()
}

private fun analyzePointMutations(transvarLocation: String, records: List<OncoActionableVariantRecord>,
                                  reference: IndexedFastaSequenceFile): List<ActionableVariantOutput> {
    val analyzer = TransvarProteinAnalyzer(transvarLocation)
    val transvarOutput = analyzer.analyze(records.map { ProteinAnnotation(it.transcript, it.alteration) })
    return records.zip(transvarOutput)
            .flatMap { (oncoRecord, transvarOutput) -> extractVariants(transvarOutput, reference).map { Pair(oncoRecord, it) } }
            .flatMap { (oncoRecord, somaticVariant) ->
                oncoRecord.drugs.map { drug ->
                    ActionableVariantOutput(oncoRecord.gene, somaticVariant,
                                            Actionability("oncoKb",
                                                          oncoRecord.cancerType,
                                                          drug,
                                                          oncoRecord.level,
                                                          oncoRecord.significance,
                                                          "predictive"))
                }
            }
}
