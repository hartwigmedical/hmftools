package com.hartwig.hmftools.knowledgebaseimporter.gene

import kotlin.math.abs
import kotlin.math.max
import kotlin.math.min

class Gene(exons: List<Exon>, private val startCodingExon: Exon, private val endCodingExon: Exon, private val seqStart: Long,
           private val seqEnd: Long, val start: Long, val end: Long, val transcript: String) {

    val chromosome = exons.first().chromosome
    private val sortedCodingExons = exons.sortedBy { it.start }.filter { it.start >= startCodingExon.start && it.start <= endCodingExon.start }
    private val exonMap: Map<Exon, Double> = annotateExons()
    val lastCodon = exonMap[sortedCodingExons.last()]!!.toInt()

    fun range(): ClosedRange<Long> = start..end

    //MIVO: returns coding ranges for the whole gene
    fun codingRanges(): List<ClosedRange<Long>> = codonCodingRanges(1, lastCodon)

    //MIVO: figures out which exon pair contains a certain codon number and returns coding ranges for that codon
    fun codonCodingRanges(codonNumber: Int): List<ClosedRange<Long>> = codonRanges(codonPositions(codonNumber))

    //MIVO: returns coding ranges between start and end codon (inclusive)
    fun codonCodingRanges(startCodon: Int, endCodon: Int): List<ClosedRange<Long>> {
        // this is different depending on forward/reverse strand
        val startCodonFirstPosition = codonPositions(startCodon)?.first
        val endCodonLastPosition = codonPositions(endCodon)?.third
        if (startCodonFirstPosition == null || endCodonLastPosition == null) return emptyList()
        return sortedCodingExons.filter { it.end >= startCodonFirstPosition && it.start <= endCodonLastPosition }
                .map { normalize(listOf(max(it.start, startCodonFirstPosition), min(it.end, endCodonLastPosition))) }
                .map { it[0]..it[1] }
                .sortedBy { it.first }
    }

    //MIVO: returns coding ranges for the specified exon number (1-based)
    fun exonCodingRanges(exonNumber: Int): List<ClosedRange<Long>> {
        val exon = sortedCodingExons.getOrNull(exonNumber - 1)
        exon ?: return emptyList()
        val normalizedPositions = normalize(listOf(codingStart(exon), codingEnd(exon)))
        return listOf(normalizedPositions[0]..normalizedPositions[1])
    }

    fun codingRangesBetween(start: Int, end: Int): List<ClosedRange<Long>> {
        return codingRanges().filter { start in it || end in it || (it.start >= start && it.endInclusive <= end) }
                .map { max(it.start, start.toLong())..min(end.toLong(), it.endInclusive) }
    }

    private fun codingStart(exon: Exon) = if (exon == startCodingExon) exon.start + seqStart - 1 else exon.start

    private fun codingEnd(exon: Exon) = if (exon == endCodingExon) exon.start + seqEnd - 1 else exon.end

    private fun codonPositions(codonNumber: Int): Triple<Long, Long, Long>? {
        if (codonNumber <= 0) return null
        val potentialExonPosition = sortedCodingExons.indexOfFirst { codonNumber <= exonMap[it]!! }
        if (potentialExonPosition == -1) return null
        val potentialExon = sortedCodingExons[potentialExonPosition]
        return if (codonOverlapsExons(codonNumber, potentialExonPosition)) {
            extractOverlappingCodon(sortedCodingExons[potentialExonPosition - 1], potentialExon)
        } else {
            extractCodonFromExon(codonNumber, potentialExon)
        }
    }

    private fun codonOverlapsExons(codonNumber: Int, potentialExonPosition: Int): Boolean {
        if (potentialExonPosition != 0) {
            val previousExon = sortedCodingExons[potentialExonPosition - 1]
            if (codonNumber == exonMap[previousExon]!!.toInt() + 1 && previousExon.endPhase > 0)
                return true
        }
        return false
    }

    //MIVO: annotate exons with the running total of codons covered so far
    //  e.g. a gene with 2 exons, each spanning 10 and 15 codons respectively, will produce:
    //      exon1 -> 10
    //      exon2 -> 25
    private fun annotateExons(): Map<Exon, Double> {
        return sortedCodingExons.fold(Pair<Long, Map<Exon, Double>>(0, mapOf())) { pair, exon ->
            val exonCodingLength = codingLength(exon)
            val runningCodingLength = pair.first + exonCodingLength
            Pair(runningCodingLength, pair.second + (exon to (runningCodingLength.toDouble() / 3)))
        }.second
    }

    //MIVO: number of coding bases contained in the exon
    private fun codingLength(exon: Exon): Long = codingEnd(exon) - codingStart(exon) + 1

    //MIVO: extract codon positions when it is contained inside an exon
    private fun extractCodonFromExon(codonNumber: Int, exon: Exon): Triple<Long, Long, Long> {
        val exonMaxCodons = exonMap[exon]!!.toInt()
        val codonsBeforeThisExon = exonMaxCodons - codingLength(exon) / 3
        val codonNumberInExon = (codonNumber - codonsBeforeThisExon - 1)
        val codonIndexInExon = codonNumberInExon * 3 + codingStart(exon) + exon.firstCodonStartOffset
        return Triple(codonIndexInExon, codonIndexInExon + 1, codonIndexInExon + 2)
    }

    //MIVO: extract codon positions when it spans 2 exons
    private fun extractOverlappingCodon(prev: Exon, next: Exon): Triple<Long, Long, Long> {
        return when {
            prev.endPhase == 1 -> Triple(prev.end, next.start, next.start + 1)
            else               -> Triple(prev.end - 1, prev.end, next.start)
        }
    }

    private fun codonRanges(triple: Triple<Long, Long, Long>?): List<ClosedRange<Long>> {
        triple ?: return emptyList()
        val positions = normalize(triple.toList())
        return when {
            positions[2] == positions[0] + 2 -> listOf(positions[0]..positions[2])
            positions[1] == positions[0] + 1 -> listOf(positions[0]..positions[1], positions[2]..positions[2])
            else                             -> listOf(positions[0]..positions[0], positions[1]..positions[2])
        }
    }

    //MIVO: convert list items to positive numbers and sort
    private fun normalize(list: List<Long>) = list.map { abs(it) }.sorted()
}
