package com.hartwig.hmftools.knowledgebaseimporter.oncoKb

import com.hartwig.hmftools.knowledgebaseimporter.output.KnownVariantOutput
import com.hartwig.hmftools.knowledgebaseimporter.transvar.TransvarProteinAnalyzer
import com.hartwig.hmftools.knowledgebaseimporter.transvar.annotations.ProteinAnnotation
import com.hartwig.hmftools.knowledgebaseimporter.transvar.extractVariants
import htsjdk.samtools.reference.IndexedFastaSequenceFile
import org.apache.commons.csv.CSVFormat
import org.apache.commons.csv.CSVParser
import java.io.File
import java.nio.charset.Charset

fun analyzeOncoKb(transvarLocation: String, fileLocation: String, reference: IndexedFastaSequenceFile): List<KnownVariantOutput> {
    val parser = CSVParser.parse(File(fileLocation), Charset.defaultCharset(), CSVFormat.TDF.withFirstRecordAsHeader().withNullString(""))
    val analyzer = TransvarProteinAnalyzer(transvarLocation)
    val records = parser.asSequence().map { OncoAnnotatedVariantRecord(it) }.toList()
    val transvarOutput = analyzer.analyze(records.map { ProteinAnnotation(it.transcript, it.alteration) })
    return records.zip(transvarOutput)
            .flatMap { (oncoRecord, transvarOutput) ->
                extractVariants(transvarOutput, reference)
                        .map { KnownVariantOutput(oncoRecord.gene, oncoRecord.transcript, oncoRecord.oncogenicity, it) }
            }
}
