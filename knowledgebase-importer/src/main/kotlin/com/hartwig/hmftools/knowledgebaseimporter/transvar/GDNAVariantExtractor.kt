package com.hartwig.hmftools.knowledgebaseimporter.transvar

import com.hartwig.hmftools.common.variant.ImmutableSomaticVariantImpl
import com.hartwig.hmftools.common.variant.SomaticVariant
import com.hartwig.hmftools.common.variant.VariantType
import htsjdk.samtools.reference.IndexedFastaSequenceFile
import org.apache.logging.log4j.LogManager
import org.apache.logging.log4j.util.Strings
import java.util.regex.Pattern

private val logger = LogManager.getLogger("GDNAVariantExtractor")

private val variantPositionsPattern = Pattern.compile("([0-9]+)(?:_([0-9]+))?.*")

private val candidateSnvString = "candidate_snv_variants="
private val candidateMnvString = "candidate_mnv_variants="
private val candidatesString = "candidates="
private val leftAlignedGDnaString = "left_align_gDNA="

private val infoDelimiter = ";"
private val variantDelimiter = ","
private val chrIdentifier = "chr"
private val chrDelimiter = ":"
private val gDnaDelimiter = "g."

fun extractVariants(transvarOutput: TransvarOutput, reference: IndexedFastaSequenceFile): List<SomaticVariant> {
    val chromosome = extractChromosome(transvarOutput.coordinates)
    return if (chromosome.isEmpty()) {
        logger.warn("could not extract chromosome from line: $transvarOutput. Skipping")
        emptyList()
    } else {
        listOfNotNull(parseOptimalCandidate(chromosome, transvarOutput, reference)) + parseInfo(chromosome, transvarOutput.info, reference)
    }
}

fun parseInfo(chromosome: String, info: String, reference: IndexedFastaSequenceFile): List<SomaticVariant> {
    val snvsGDna = info.substringAfter(candidateSnvString,
                                       "").substringBefore(infoDelimiter).split(variantDelimiter).filterNot { it.isEmpty() }
    val mnvsGDna = info.substringAfter(candidateMnvString,
                                       "").substringBefore(infoDelimiter).split(variantDelimiter).filterNot { it.isEmpty() }
    val indelsGDna = info.substringAfter(candidatesString,
                                         "").substringBefore(infoDelimiter).split(variantDelimiter).filterNot { it.isEmpty() }.map {
        it.split("/")[2]
    }
    val mergedVariants = snvsGDna + mnvsGDna + indelsGDna
    return extract(chromosome, mergedVariants, reference)
}

fun parseOptimalCandidate(chromosome: String, transvarOutput: TransvarOutput, reference: IndexedFastaSequenceFile): SomaticVariant? {
    val variantGDna = if (transvarOutput.info.contains(leftAlignedGDnaString)) {
        transvarOutput.info.substringAfter(leftAlignedGDnaString).substringAfter(gDnaDelimiter, "").substringBefore(";")
    } else {
        transvarOutput.coordinates.substringAfter(gDnaDelimiter, "").substringBefore("/")
    }
    return extractVariant(chromosome, variantGDna, reference)
}

fun extract(chromosome: String, gDnaVariants: List<String>, reference: IndexedFastaSequenceFile): List<SomaticVariant> {
    return gDnaVariants.map { it.substringAfter(gDnaDelimiter) }.mapNotNull {
        extractVariant(chromosome, it, reference)
    }
}

fun extractVariant(chromosome: String, variantGDna: String, reference: IndexedFastaSequenceFile): SomaticVariant? {
    return try {
        when {
            variantGDna.contains(">")                                  -> extractSnv(chromosome, variantGDna)
            variantGDna.contains("del") && variantGDna.contains("ins") -> extractMnv(chromosome, variantGDna, reference)
            variantGDna.contains("ins")                                -> extractInsert(chromosome, variantGDna, reference)
            variantGDna.contains("dup")                                -> extractDup(chromosome, variantGDna, reference)
            variantGDna.contains("del")                                -> extractDelete(chromosome, variantGDna, reference)
            else                                                       -> {
                logger.warn("variant $chromosome: $variantGDna could not be mapped to any known type")
                null
            }
        }
    } catch (t: Throwable) {
        logger.warn("Could not create variant from $chromosome: $variantGDna; error: $t")
        null
    }
}

fun extractSnv(chromosome: String, variantGDna: String): SomaticVariant {
    val (start, _) = extractPositions(variantGDna)
    val ref = variantGDna.substringBefore(">").last().toString()
    val alt = variantGDna.substringAfter(">").first().toString()
    return somaticVariant(chromosome, start, ref, alt)
}

fun extractChromosome(coordinates: String): String {
    return coordinates.substringAfter(chrIdentifier, "").substringBefore(chrDelimiter)
}

private fun extractPositions(variant: String): Pair<Long, Long?> {
    val matcher = variantPositionsPattern.matcher(variant.substringAfter(gDnaDelimiter))
    matcher.find()
    return Pair(matcher.group(1).toLong(), matcher.group(2)?.toLongOrNull())
}

fun extractMnv(chromosome: String, variantGDna: String, reference: IndexedFastaSequenceFile): SomaticVariant {
    val (start, end) = extractPositions(variantGDna)
    return if (variantGDna.contains("delins")) {
        val endPosition = end ?: start
        val alt = variantGDna.substringAfter("delins")
        val ref = reference.getSubsequenceAt(chromosome, start, endPosition).baseString
        somaticVariant(chromosome, start, ref, alt)
    } else {
        val ref = variantGDna.substringAfter("del").substringBefore("ins")
        val alt = variantGDna.substringAfter("ins")
        return somaticVariant(chromosome, start, ref, alt)
    }
}

fun extractDup(chromosome: String, variantGDna: String, reference: IndexedFastaSequenceFile): SomaticVariant {
    val (start, _) = extractPositions(variantGDna)
    val position = start - 1
    val ref = reference.getSubsequenceAt(chromosome, position, position).baseString
    val insertedBases = variantGDna.substringAfter("dup")
    val alt = ref + insertedBases
    return somaticVariant(chromosome, position, ref, alt)
}

fun extractInsert(chromosome: String, variantGDna: String, reference: IndexedFastaSequenceFile): SomaticVariant? {
    val (start, _) = extractPositions(variantGDna)
    val ref = reference.getSubsequenceAt(chromosome, start, start).baseString
    val insertedBases = variantGDna.substringAfter("ins")
    val alt = ref + insertedBases
    return if (alt.contains("N")) {
        null
    } else {
        somaticVariant(chromosome, start, ref, alt)
    }
}

fun extractDelete(chromosome: String, variantGDna: String, reference: IndexedFastaSequenceFile): SomaticVariant? {
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
            logger.warn("Skipping deletion of more than 20 bases for variant $chromosome: $variantGDna")
            return null
        }
        val ref = reference.getSubsequenceAt(chromosome, position, position + deletedBasesCount).baseString
        val alt = ref.first().toString()
        Pair(ref, alt)
    }
    return somaticVariant(chromosome, position, ref, alt)
}

fun somaticVariant(chromosome: String, position: Long, ref: String, alt: String): SomaticVariant {
    //type, filter, gene, genesEffected, worstEffect, worstEffectTranscript, worstCodingEffect, hotspot, mappability, totalReadCount, alleleReadCount
    return ImmutableSomaticVariantImpl.builder().chromosome(chromosome).position(position).ref(ref).alt(alt)
            .type(VariantType.UNDEFINED)
            .filter(Strings.EMPTY)
            .gene(Strings.EMPTY)
            .genesEffected(0)
            .worstEffect(Strings.EMPTY)
            .worstEffectTranscript(Strings.EMPTY)
            .worstCodingEffect(Strings.EMPTY)
            .totalReadCount(0)
            .alleleReadCount(0)
            .hotspot(false)
            .mappability(0.0).build()
}
