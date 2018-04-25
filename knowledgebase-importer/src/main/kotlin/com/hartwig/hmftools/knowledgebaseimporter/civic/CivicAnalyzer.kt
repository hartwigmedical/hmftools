package com.hartwig.hmftools.knowledgebaseimporter.civic

import com.hartwig.hmftools.common.variant.SomaticVariant
import com.hartwig.hmftools.knowledgebaseimporter.output.Actionability
import com.hartwig.hmftools.knowledgebaseimporter.output.ActionableCNVOutput
import com.hartwig.hmftools.knowledgebaseimporter.output.ActionableVariantOutput
import com.hartwig.hmftools.knowledgebaseimporter.output.KnownVariantOutput
import com.hartwig.hmftools.knowledgebaseimporter.transvar.TransvarCdnaAnalyzer
import com.hartwig.hmftools.knowledgebaseimporter.transvar.annotations.CDnaAnnotation
import com.hartwig.hmftools.knowledgebaseimporter.transvar.extractVariants
import com.hartwig.hmftools.knowledgebaseimporter.transvar.somaticVariant
import htsjdk.samtools.reference.IndexedFastaSequenceFile

private const val SOURCE = "civic"

fun readCivicVariants(transvarLocation: String, variantFileLocation: String, evidenceFileLocation: String,
                      reference: IndexedFastaSequenceFile): List<Pair<CivicRecord, SomaticVariant>> {
    val records = preProcessCivicVariants(variantFileLocation, evidenceFileLocation)
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
                (listOfNotNull(civicVariant) + inferredVariants.filterNot { it == civicVariant }).map { Pair(civicRecord, it) }
            }
}

fun analyzeKnownCivicVariants(transvarLocation: String, variantFileLocation: String, evidenceFileLocation: String,
                              reference: IndexedFastaSequenceFile): List<KnownVariantOutput> {
    val civicVariants = readCivicVariants(transvarLocation, variantFileLocation, evidenceFileLocation, reference)
    return civicVariants.map { (civicRecord, somaticVariant) ->
        KnownVariantOutput(civicRecord.gene, civicRecord.transcript, additionalInfo(civicRecord), somaticVariant)
    }
}

fun analyzeCivicActionable(transvarLocation: String, variantFileLocation: String, evidenceFileLocation: String,
                           reference: IndexedFastaSequenceFile): List<ActionableVariantOutput> {
    val civicVariants = readCivicVariants(transvarLocation, variantFileLocation, evidenceFileLocation, reference)
    return civicVariants.flatMap { (record, somaticVariant) ->
        record.evidence.filter { it.direction == "Supports" }.flatMap { evidence ->
            evidence.drugs.map { drug ->
                ActionableVariantOutput(record.gene,
                                        somaticVariant,
                                        Actionability(SOURCE,
                                                      evidence.cancerType,
                                                      drug,
                                                      evidence.level,
                                                      evidence.significance,
                                                      evidence.type))
            }
        }
    }
}

private fun additionalInfo(civicRecord: CivicRecord): String {
    val highestEvidenceLevel = civicRecord.evidence.map { it.level }.sorted().firstOrNull() ?: "N"
    return (highestEvidenceLevel == "A" || highestEvidenceLevel == "B" || highestEvidenceLevel == "C").toString()
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

private fun extractAmpOrDel(variantType: String): String = if (variantType == "AMPLIFICATION") "Amplification" else "Deletion"

fun analyzeCivicAmpsAndDels(variantFileLocation: String, evidenceFileLocation: String): List<ActionableCNVOutput> {
    val records = preProcessCivicRecords(variantFileLocation, evidenceFileLocation)
    return records.filter { it.variant == "AMPLIFICATION" || it.variant == "DELETION" || it.variant == "LOH" }.flatMap { record ->
        record.evidence.flatMap { evidence ->
            evidence.drugs.map { drug ->
                ActionableCNVOutput(record.gene, extractAmpOrDel(record.variant),
                                    Actionability(SOURCE, evidence.cancerType, drug, evidence.level, evidence.significance, evidence.type))
            }
        }
    }
}
