package com.hartwig.hmftools.lilac

import com.hartwig.hmftools.lilac.amino.AminoAcidFragment
import com.hartwig.hmftools.lilac.nuc.NucleotideFragment
import java.io.File
import java.util.*
import kotlin.math.min

class SequenceCount(private val minCount: Int, val length: Int) {

    private val count = Array(length) { mutableMapOf<String, Int>() }

    companion object {
        fun nucleotides(minCount: Int, fragments: List<NucleotideFragment>): SequenceCount {
            val length = fragments.map { it.nucleotideLoci().max() ?: -1 }.max()!! + 1
            val result = SequenceCount(minCount, length)

            for (fragment in fragments) {
                for (index in fragment.nucleotideLoci()) {
                    val nucleotide = fragment.nucleotide(index)
                    result.increment(index, nucleotide)
                }
            }
            return result
        }

        fun aminoAcids(minCount: Int, aminoAcidFragments: List<AminoAcidFragment>): SequenceCount {
            val length = aminoAcidFragments.map { it.aminoAcidIndices().max() ?: -1 }.max()!! + 1
            val result = SequenceCount(minCount, length)

            for (fragment in aminoAcidFragments) {
                for (index in fragment.aminoAcidIndices()) {
                    val aminoAcid = fragment.aminoAcid(index)
                        result.increment(index, aminoAcid)
                }
            }
            return result
        }
    }

    operator fun get(locus: Int): Map<String, Int>  {
        return count[locus]
    }

    private fun increment(index: Int, aminoAcid: String) {
        try {
            count[index].compute(aminoAcid) { _, u -> (u ?: 0) + 1 }
        } catch (e: Exception) {
            throw e
        }
    }

    fun heterozygousLoci(): List<Int> {
        return count.indices.filter { isHeterozygous(it) }
    }

    fun homozygousIndices(): List<Int> {
        return count.indices.filter { isHomozygous(it) }
    }

    private fun isHomozygous(index: Int): Boolean {
        val mapAtIndex = count[index].filter { it.value >= minCount }
        return mapAtIndex.size == 1
    }


    private fun isHeterozygous(index: Int): Boolean {
        val mapAtIndex = count[index].filter { it.value >= minCount }
        return mapAtIndex.size > 1
    }

    fun sequenceAt(index: Int): Collection<String> {
        if (index >= count.size) {
            return Collections.emptySet()
        }

        val result = mutableSetOf<String>()
        val indexMap = count[index]
        for ((aa, count) in indexMap) {
            if (count >= minCount) {
                result.add(aa)
            }
        }
        return result
    }


    fun depth(index: Int): Int {
        val indexMap = count[index]
        return indexMap.values.sum()
    }


    fun writeVertically(fileName: String) {
        val file = File(fileName)
        file.writeText("")
        for (i in 0 until length) {
            val lineBuilder = StringJoiner("\t").add(i.toString())
            val baseCountList = count[i].map { (k, v) -> Pair(k, v) }.sortedBy { it.second }.reversed()
            for (j in 0..min(5, baseCountList.size - 1)) {
                val (base, count) = baseCountList[j]
                lineBuilder.add(base.toString()).add(count.toString())
            }

            file.appendText(lineBuilder.toString() + "\n")
        }
    }

}