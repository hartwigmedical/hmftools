package com.hartwig.hmftools.lilac.coverage

import com.hartwig.hmftools.common.utils.SuffixTree
import com.hartwig.hmftools.lilac.dna.aminoAcids
import com.hartwig.hmftools.lilac.dna.dnaReverseComplement
import htsjdk.samtools.SAMRecord
import java.util.function.Consumer

class BamCoverage(val minMatch: Int, val proteins: Set<ProteinCoverage>) : Consumer<SAMRecord> {

    override fun accept(record: SAMRecord) {
        // Do some quality trimming???
        processDna(record.readString)
    }

    fun processDna(dna: String) {
        val dnaReversed = dna.dnaReverseComplement()

        val codonString0 = dna.aminoAcids()
        val codonString1 = dna.substring(1).aminoAcids()
        val codonString2 = dna.substring(2).aminoAcids()

        val codonString3 = dnaReversed.aminoAcids()
        val codonString4 = dnaReversed.substring(1).aminoAcids()
        val codonString5 = dnaReversed.substring(2).aminoAcids()

        val searchTrees = listOf(
                SuffixTree(codonString0),
                SuffixTree(codonString1),
                SuffixTree(codonString2),
                SuffixTree(codonString3),
                SuffixTree(codonString4),
                SuffixTree(codonString5)
        )

        val reverseSearchTrees = listOf(
                SuffixTree(codonString0.reversed()),
                SuffixTree(codonString1.reversed()),
                SuffixTree(codonString2.reversed()),
                SuffixTree(codonString3.reversed()),
                SuffixTree(codonString4.reversed()),
                SuffixTree(codonString5.reversed()))

        val (endMatchSize, endMatchExons) = largestMatch(searchTrees)
        val (startMatchSize, startMatchExons) = largestMatch(reverseSearchTrees, true)
        if (endMatchSize >= minMatch || startMatchSize >= minMatch) {
            when {
                endMatchSize > startMatchSize -> endMatchExons.forEach { it.addCoverageFromStart(endMatchSize) }
                endMatchSize < startMatchSize -> startMatchExons.forEach { it.addCoverageFromEnd(startMatchSize) }
                else -> {
                    endMatchExons
                            .forEach { it.addCoverageFromStart(endMatchSize) }
                    startMatchExons
                            .filter { !endMatchExons.contains(it) } // Don't double up!
                            .forEach { it.addCoverageFromEnd(startMatchSize) }
                }
            }
        }
    }


    private fun largestMatch(searchTrees: Collection<SuffixTree>, reverse: Boolean = false): Pair<Int, List<ExonCoverage>> {
        var largestMatch = 0
        val largest = mutableListOf<ExonCoverage>()
        for (protein in proteins) {
            for ((exon, coverage) in protein.map.entries) {
                val exonToMatch = if (reverse) exon.reversed() else exon

                for (tree in searchTrees) {
                    val matchingBases = if (tree.contains(exonToMatch)) exon.length else tree.endsWith(exonToMatch)
                    if (matchingBases > largestMatch) {
                        largest.clear()
                        largestMatch = matchingBases
                    }

                    if (matchingBases == largestMatch) {
                        largest.add(coverage)
                    }
                }
            }
        }

        return Pair(largestMatch, largest)
    }

}