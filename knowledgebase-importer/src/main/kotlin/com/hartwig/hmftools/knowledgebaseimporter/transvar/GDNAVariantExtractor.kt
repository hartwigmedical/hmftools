package com.hartwig.hmftools.knowledgebaseimporter.transvar

import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.KnowledgebaseRecord
import com.hartwig.hmftools.knowledgebaseimporter.output.ActionableEvent
import com.hartwig.hmftools.knowledgebaseimporter.output.GenomicRangeEvent
import com.hartwig.hmftools.knowledgebaseimporter.output.SomaticVariantEvent
import htsjdk.samtools.reference.IndexedFastaSequenceFile
import org.apache.logging.log4j.LogManager
import java.util.regex.Pattern

private val logger = LogManager.getLogger("GDNAVariantExtractor")

private val variantPositionsPattern = Pattern.compile("([0-9]+)(?:_([0-9]+))?.*")

private const val candidateSnvString = "candidate_snv_variants="
private const val candidateMnvString = "candidate_mnv_variants="
private const val candidatesString = "candidates="
private const val leftAlignedGDnaString = "left_align_gDNA="

private const val infoDelimiter = ";"
private const val variantDelimiter = ","
private const val chrIdentifier = "chr"
private const val chrDelimiter = ":"
private const val gDnaDelimiter = "g."

fun extractVariants(record: KnowledgebaseRecord, transvarOutput: TransvarOutput, reference: IndexedFastaSequenceFile):
        List<ActionableEvent> {
    return if (transvarOutput.info == "no_valid_transcript_found") {
        logger.warn("Transvar could not resolve genomic coordinates for gene ${record.gene} with events ${record.events}" )
        emptyList()
    } else {
        val chromosome = extractChromosome(transvarOutput.coordinates)
        return if (chromosome.isEmpty()) {
            logger.warn("Could not extract chromosome for $record from transvar output: $transvarOutput. Skipping")
            emptyList()
        } else {
            listOfNotNull(parseOptimalCandidate(record, chromosome, transvarOutput, reference)) +
                    parseInfo(record, chromosome, transvarOutput.info, reference)
        }
    }
}

private fun parseInfo(record: KnowledgebaseRecord, chromosome: String, info: String,
                      reference: IndexedFastaSequenceFile): List<ActionableEvent> {
    val snvsGDna = infoVariantsAfter(info, candidateSnvString)
    val mnvsGDna = infoVariantsAfter(info, candidateMnvString)
    val indelsGDna = infoVariantsAfter(info, candidatesString).map { it.split("/")[2] }
    val mergedVariants = snvsGDna + mnvsGDna + indelsGDna
    return extractVariants(record, chromosome, mergedVariants, reference)
}

private fun infoVariantsAfter(info: String, startDelimiter: String): List<String> {
    return info.substringAfter(startDelimiter, "").substringBefore(infoDelimiter).split(variantDelimiter).filterNot { it.isEmpty() }
}

private fun parseOptimalCandidate(record: KnowledgebaseRecord, chromosome: String, transvarOutput: TransvarOutput,
                                  reference: IndexedFastaSequenceFile): ActionableEvent? {
    val variantGDna = if (transvarOutput.info.contains(leftAlignedGDnaString)) {
        transvarOutput.info.substringAfter(leftAlignedGDnaString).substringAfter(gDnaDelimiter, "").substringBefore(";")
    } else {
        transvarOutput.coordinates.substringAfter(gDnaDelimiter, "").substringBefore("/")
    }
    return extractVariant(record, chromosome, variantGDna, reference)
}

private fun extractVariants(record: KnowledgebaseRecord, chromosome: String, gDnaVariants: List<String>,
                            reference: IndexedFastaSequenceFile): List<ActionableEvent> {
    return gDnaVariants.map { it.substringAfter(gDnaDelimiter) }.mapNotNull {
        extractVariant(record, chromosome, it, reference)
    }
}

fun extractVariant(record: KnowledgebaseRecord, gDnaVariant: String, reference: IndexedFastaSequenceFile): ActionableEvent? {
    val chromosome = extractChromosome(gDnaVariant)
    return extractVariant(record, chromosome, gDnaVariant.substringAfter(gDnaDelimiter), reference)
}

private fun extractVariant(record: KnowledgebaseRecord, chromosome: String, variantGDna: String,
                           reference: IndexedFastaSequenceFile): ActionableEvent? {
    return try {
        val gene = record.gene
        val transcript = record.transcript
        when {
            variantGDna.contains(">") -> extractSnv(gene, chromosome, variantGDna)
            variantGDna.contains("del") && variantGDna.contains("ins") -> extractMnv(gene, chromosome, variantGDna, reference)
            variantGDna.contains("ins") -> extractInsert(gene, chromosome, variantGDna, reference)
            variantGDna.contains("dup") -> extractDup(gene, chromosome, variantGDna, reference)
            variantGDna.contains("del") -> extractDelete(gene, chromosome, variantGDna, reference)
            variantGDna.matches("[0-9]+_[0-9]+".toRegex()) -> extractRange(gene, transcript, chromosome, variantGDna)
            else -> {
                logger.warn("Gene ${record.gene} with events ${record.events} " +
                        "could not be mapped to any known type based on variant gDNA $variantGDna")
                null
            }
        }
    } catch (t: Throwable) {
        logger.warn("Could not create variant on ${record.gene} on position $chromosome:$variantGDna; error: $t")
        null
    }
}

fun extractChromosome(coordinates: String): String {
    return coordinates.substringAfter(chrIdentifier, "").substringBefore(chrDelimiter)
}

private fun extractPositions(variant: String): Pair<Long, Long?> {
    val matcher = variantPositionsPattern.matcher(variant.substringAfter(gDnaDelimiter))
    matcher.find()
    return Pair(matcher.group(1).toLong(), matcher.group(2)?.toLongOrNull())
}

// MIVO: extract SNV from gDna variant of the form: 133738357T>C
private fun extractSnv(gene: String, chromosome: String, variantGDna: String): SomaticVariantEvent {
    val (start, _) = extractPositions(variantGDna)
    val ref = variantGDna.substringBefore(">").last().toString()
    val alt = variantGDna.substringAfter(">").first().toString()
    return SomaticVariantEvent(gene, chromosome, start.toString(), ref, alt)
}

// MIVO: extract MNV from gDna variant of the form: 105239404_105239405delinsGC or 133748289_133748290delTCinsGT
private fun extractMnv(gene: String, chromosome: String, variantGDna: String, reference: IndexedFastaSequenceFile): SomaticVariantEvent {
    val (start, end) = extractPositions(variantGDna)
    return if (variantGDna.contains("delins")) {
        val endPosition = end ?: start
        val alt = variantGDna.substringAfter("delins")
        val ref = reference.getSubsequenceAt(chromosome, start, endPosition).baseString
        correct(SomaticVariantEvent(gene, chromosome, start.toString(), ref, alt))
    } else {
        val ref = variantGDna.substringAfter("del").substringBefore("ins")
        val alt = variantGDna.substringAfter("ins")
        return correct(SomaticVariantEvent(gene, chromosome, start.toString(), ref, alt))
    }
}

// MIVO: extract Insert from gDna variant of the form: 41201160dupA
private fun extractDup(gene: String, chromosome: String, variantGDna: String, reference: IndexedFastaSequenceFile): SomaticVariantEvent {
    val (start, _) = extractPositions(variantGDna)
    val position = start - 1
    val ref = reference.getSubsequenceAt(chromosome, position, position).baseString
    val insertedBases = variantGDna.substringAfter("dup")
    val alt = ref + insertedBases
    return correct(SomaticVariantEvent(gene, chromosome, position.toString(), ref, alt))
}

// MIVO: extract Insert from gDna variant of the form: 32930598_32930599insC
private fun extractInsert(gene: String, chromosome: String, variantGDna: String,
                          reference: IndexedFastaSequenceFile): SomaticVariantEvent? {
    val (start, _) = extractPositions(variantGDna)
    val ref = reference.getSubsequenceAt(chromosome, start, start).baseString
    val insertedBases = variantGDna.substringAfter("ins")
    val alt = ref + insertedBases
    return if (alt.contains("N")) {
        null
    } else {
        correct(SomaticVariantEvent(gene, chromosome, start.toString(), ref, alt))
    }
}

// MIVO: extract Delete from gDna variant of the form: 55152094_55152105delCATCATGCATGA or 55152095_55152106del12
private fun extractDelete(gene: String, chromosome: String, variantGDna: String,
                          reference: IndexedFastaSequenceFile): SomaticVariantEvent? {
    val (start, _) = extractPositions(variantGDna)
    val position = start - 1
    val deletedBases = variantGDna.substringAfter("del")
    val deletedBasesCount = deletedBases.toIntOrNull()
    val (ref, alt) = if (deletedBasesCount == null) {
        val alt = reference.getSubsequenceAt(chromosome, position, position).baseString
        val ref = alt + deletedBases
        Pair(ref, alt)
    } else {
        if (deletedBasesCount > 20) {
            logger.info("Skipping deletion of more than 20 bases for variant on $gene on position $chromosome:$variantGDna")
            return null
        }
        val ref = reference.getSubsequenceAt(chromosome, position, position + deletedBasesCount).baseString
        val alt = ref.first().toString()
        Pair(ref, alt)
    }
    return correct(SomaticVariantEvent(gene, chromosome, position.toString(), ref, alt))
}

private fun extractRange(gene: String, transcript: String, chromosome: String, variantGDna: String): GenomicRangeEvent {
    val (start, end) = extractPositions(variantGDna)
    return GenomicRangeEvent(gene, transcript, chromosome, start.toString(), end.toString(), transcript)
}

// MIVO: correct position, ref and alt of variants for which ref and alt have a common prefix.
// e.g. chr3:g.178927980_178927991delinsTGG will produce: 3 178927980: TGTCCATTGGCA -> TGG, with both ref and alt sharing TG as prefix
//      this function will drop the first T and correct the variant to: 3 178927981: GTCCATTGGCA -> GG
private fun correct(variant: SomaticVariantEvent): SomaticVariantEvent {
    val equalBaseCount = variant.ref.zip(variant.alt).takeWhile { it.first == it.second }.count()
    return if (equalBaseCount > 1) {
        val basesToSkip = equalBaseCount - 1
        val position = variant.position.toInt() + basesToSkip
        variant.copy(position = position.toString(), ref = variant.ref.substring(basesToSkip), alt = variant.alt.substring(basesToSkip))
    } else {
        variant
    }
}
