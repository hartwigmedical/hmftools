package com.hartwig.hmftools.knowledgebaseimporter.dao

import kotlin.math.abs

class Gene(exons: List<Exon>, private val seqStart: Long, private val seqEnd: Long) {

    val sortedCodingExons = exons.filterNot { it.isNonCoding }.sortedBy { it.start }
    private val exonMap: Map<Exon, Double> = annotateExons()

    //MIVO: figures out which exon pair contains a certain codon number and extracts codon positions
    fun codonPositions(codonNumber: Int): Triple<Long, Long, Long>? {
        if (codonNumber <= 0) return null
        return sortedCodingExons.zipWithNext().find { (prev, next) ->
            exonMap[prev]!! >= codonNumber || exonMap[next]!! >= codonNumber
        }?.run { extractCodonPositions(codonNumber, first, second) }
    }

    //MIVO: annotate exons with the running total of codons covered so far
    //  e.g. a gene with 2 exons, each spanning 10 and 15 codons respectively, will produce:
    //      exon1 -> 10
    //      exon2 -> 25
    fun annotateExons(): Map<Exon, Double> {
        return sortedCodingExons.fold(Pair<Long, Map<Exon, Double>>(0, mapOf())) { pair, exon ->
            val exonCodingLength = codingLength(exon)
            val runningCodingLength = pair.first + exonCodingLength
            Pair(runningCodingLength, pair.second + (exon to (runningCodingLength.toDouble() / 3)))
        }.second
    }

    //MIVO: number of coding bases contained in the exon
    private fun codingLength(exon: Exon): Long {
        return when {
            exon.isFirst -> exon.length - seqStart + 1
            exon.isLast  -> seqEnd
            else         -> exon.length
        }
    }

    private fun extractCodonPositions(codonNumber: Int, prevExon: Exon, nextExon: Exon): Triple<Long, Long, Long> {
        val prevMaxCodons = exonMap[prevExon]!!.toInt()
        return when {
            codonNumber <= prevMaxCodons                                -> extractCodonFromExon(codonNumber, prevExon)
            codonNumber == (prevMaxCodons + 1) && prevExon.endPhase > 0 -> extractOverlappingCodon(prevExon, nextExon)
            else                                                        -> extractCodonFromExon(codonNumber, nextExon)
        }
    }

    //MIVO: extract codon positions when it is contained inside an exon
    private fun extractCodonFromExon(codonNumber: Int, exon: Exon): Triple<Long, Long, Long> {
        val exonMaxCodons = exonMap[exon]!!.toInt()
        val codonsBeforeThisExon = exonMaxCodons - codingLength(exon) / 3
        val codonNumberInExon = (codonNumber - codonsBeforeThisExon - 1)
        val seqStartOffset = if (exon.isFirst) seqStart - 1 else 0
        val codonIndexInExon = codonNumberInExon * 3 + exon.start + exon.firstCodonStartOffset + seqStartOffset
        return codonPositions(codonIndexInExon, codonIndexInExon + 1, codonIndexInExon + 2)
    }

    //MIVO: extract codon positions when it spans 2 exons
    private fun extractOverlappingCodon(prev: Exon, next: Exon): Triple<Long, Long, Long> {
        return when {
            prev.endPhase == 1 -> codonPositions(prev.end, next.start, next.start + 1)
            else               -> codonPositions(prev.end - 1, prev.end, next.start)
        }
    }

    //MIVO: normalize codon positions (reverse strand positions are converted to positive numbers and sorted)
    private fun codonPositions(first: Long, second: Long, third: Long): Triple<Long, Long, Long> {
        val positions = listOf(first, second, third).map { abs(it) }.sorted()
        return Triple(positions[0], positions[1], positions[2])
    }
}
