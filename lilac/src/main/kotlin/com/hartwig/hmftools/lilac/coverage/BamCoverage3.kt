package com.hartwig.hmftools.lilac.coverage

import com.hartwig.hmftools.common.utils.SuffixTree
import com.hartwig.hmftools.lilac.dna.aminoAcids
import com.hartwig.hmftools.lilac.dna.dnaReverseComplement
import htsjdk.samtools.SAMRecord
import java.util.function.Consumer

class BamCoverage3(val minOverlap: Int, coverage: Collection<ExonCoverage>) : Consumer<SAMRecord> {

    companion object {
        const val TRIM = 10
    }

    private val forwardTrees: Map<ExonCoverage, SuffixTree> = coverage.map { Pair(it, SuffixTree(it.exonSequence)) }.toMap()
    private val reverseTrees: Map<ExonCoverage, SuffixTree> = coverage.map { Pair(it, SuffixTree(it.exonSequence.reversed())) }.toMap()

    override fun accept(record: SAMRecord) {
        val readString = record.readString
        processDna(readString.substring(TRIM, readString.length - TRIM))
    }

    fun processDna(dna: String) {

        val dnaReversed = dna.dnaReverseComplement()

        val codonString0 = dna.aminoAcids()
        val codonString1 = dna.substring(1).aminoAcids()
        val codonString2 = dna.substring(2).aminoAcids()

        val codonString3 = dnaReversed.aminoAcids()
        val codonString4 = dnaReversed.substring(1).aminoAcids()
        val codonString5 = dnaReversed.substring(2).aminoAcids()

        val allCodonStrings = listOf(codonString0, codonString1, codonString2, codonString3, codonString4, codonString5)


        val endsWith = { tree: SuffixTree ->
            intArrayOf(
                    tree.endsWith(codonString0),
                    tree.endsWith(codonString1),
                    tree.endsWith(codonString2),
                    tree.endsWith(codonString3),
                    tree.endsWith(codonString4),
                    tree.endsWith(codonString5)).max()!!
        }

        val startsWith = { tree: SuffixTree ->
            intArrayOf(
                    tree.endsWith(codonString0.reversed()),
                    tree.endsWith(codonString1.reversed()),
                    tree.endsWith(codonString2.reversed()),
                    tree.endsWith(codonString3.reversed()),
                    tree.endsWith(codonString4.reversed()),
                    tree.endsWith(codonString5.reversed())).max()!!
        }

        // Check contains
        var foundExactMatch = false
        for ((exonCoverage, tree) in forwardTrees) {

            for (codonString in allCodonStrings) {
                val index = tree.anyIndexOf(codonString)
                if (index != -1) {
                    foundExactMatch = true
                    exonCoverage.addCoverage(index, codonString.length)
                    break
                }
            }
        }

        if (foundExactMatch) {
            return
        }

        val (longestEndLength, longestEnds) = overlap(forwardTrees, endsWith)
        val (longestStartLength, longestStarts) = overlap(reverseTrees, startsWith)

        if (longestEndLength >= minOverlap && longestEndLength >= longestStartLength) {
            longestEnds.forEach { it.addCoverageFromEnd(longestEndLength) }
        }

        if (longestStartLength >= minOverlap && longestStartLength >= longestEndLength) {
            longestStarts.forEach { it.addCoverageFromStart(longestStartLength) }
        }
    }


    private fun overlap(map: Map<ExonCoverage, SuffixTree>, coverageFunction: (SuffixTree) -> Int): Pair<Int, List<ExonCoverage>> {
        var longest = 0
        val result = mutableListOf<ExonCoverage>()
        for ((exonCoverage, tree) in map) {
            val coverage = coverageFunction(tree)
            if (coverage > longest && coverage >= minOverlap) {
                result.clear()
                longest = coverage
            }

            if (coverage == longest) {
                result.add(exonCoverage)
            }
        }

        return Pair(longest, result)
    }
}