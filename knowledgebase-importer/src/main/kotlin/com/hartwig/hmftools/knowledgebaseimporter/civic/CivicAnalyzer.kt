package com.hartwig.hmftools.knowledgebaseimporter.civic

import com.hartwig.hmftools.common.variant.SomaticVariant
import com.hartwig.hmftools.knowledgebaseimporter.output.KnownVariantOutput
import com.hartwig.hmftools.knowledgebaseimporter.transvar.TransvarCdnaAnalyzer
import com.hartwig.hmftools.knowledgebaseimporter.transvar.annotations.CDnaAnnotation
import com.hartwig.hmftools.knowledgebaseimporter.transvar.extractVariants
import com.hartwig.hmftools.knowledgebaseimporter.transvar.somaticVariant
import htsjdk.samtools.reference.IndexedFastaSequenceFile

fun analyzeCivic(transvarLocation: String, variantFileLocation: String, evidenceFileLocation: String,
                 reference: IndexedFastaSequenceFile): List<KnownVariantOutput> {
    val records = preProcessCivic(variantFileLocation, evidenceFileLocation)
    val analyzer = TransvarCdnaAnalyzer(transvarLocation)
    val cdnaRecords = records.map { record ->
        val hgvsParts = record.hgvs.split(":")
        if (hgvsParts.size > 1) {
            CDnaAnnotation(hgvsParts[0], hgvsParts[1])
        } else {
            CDnaAnnotation("na", "na")
        }
    }
    val transvarOutput = analyzer.analyze(cdnaRecords)
    return records.zip(transvarOutput)
            .flatMap { (civicRecord, transvarOutput) ->
                val civicVariant = annotateCivicVariant(civicRecord, reference)
                val inferredVariants = extractVariants(transvarOutput, reference)
                (listOfNotNull(civicVariant) + inferredVariants.filterNot { it == civicVariant }).map {
                    KnownVariantOutput(civicRecord.gene, civicRecord.transcript, additionalInfo(civicRecord), it)
                }
            }
}

private fun additionalInfo(civicRecord: CivicRecord): String {
    val evidenceLevel = civicRecord.evidenceLevel
    return (evidenceLevel == "A" || evidenceLevel == "B" || evidenceLevel == "C").toString()
}

private fun annotateCivicVariant(civicRecord: CivicRecord, reference: IndexedFastaSequenceFile): SomaticVariant? {
    return if (!civicRecord.chromosome.isEmpty() && !civicRecord.start.isEmpty() && (!civicRecord.ref.isEmpty() || !civicRecord.alt.isEmpty())) {
        val position = civicRecord.start.toLong()
        when {
            civicRecord.ref.isEmpty() -> {
                val base = reference.getSubsequenceAt(civicRecord.chromosome, position, position).baseString
                somaticVariant(civicRecord.chromosome, position, base, base + civicRecord.alt)
            }
            civicRecord.alt.isEmpty() -> {
                val base = reference.getSubsequenceAt(civicRecord.chromosome, position - 1, position - 1).baseString
                somaticVariant(civicRecord.chromosome, position - 1, base + civicRecord.ref, base)
            }
            else                      -> somaticVariant(civicRecord.chromosome, position, civicRecord.ref, civicRecord.alt)
        }
    } else {
        null
    }
}
