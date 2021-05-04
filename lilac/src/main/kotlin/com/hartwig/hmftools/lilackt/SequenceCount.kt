package com.hartwig.hmftools.lilackt

import com.hartwig.hmftools.lilackt.amino.AminoAcidFragment
import com.hartwig.hmftools.lilackt.nuc.NucleotideFragment
import java.io.File
import java.util.*
import kotlin.math.min

class SequenceCount(private val minCount: Int, private val count: Array<Map<String, Int>>) {
    val length: Int = count.size

    companion object {

        fun nucleotides(minCount: Int, fragments: List<NucleotideFragment>): SequenceCount {
            val length = fragments.map { it.nucleotideLoci().max() ?: -1 }.max()!! + 1
            val count = Array(length) { mutableMapOf<String, Int>() }

            for (fragment in fragments) {
                for (index in fragment.nucleotideLoci()) {
                    val nucleotide = fragment.nucleotide(index)
                    count.increment(index, nucleotide)
                }
            }
            return SequenceCount(minCount, Array(length) { count[it] })
        }

        fun aminoAcids(minCount: Int, fragments: List<AminoAcidFragment>): SequenceCount {
            val length = fragments.map { it.aminoAcidLoci().max() ?: -1 }.max()!! + 1
            val count = Array(length) { mutableMapOf<String, Int>() }

            for (fragment in fragments) {
                for (index in fragment.aminoAcidLoci()) {
                    val aminoAcid = fragment.aminoAcid(index)
                    count.increment(index, aminoAcid)
                }
            }
            return SequenceCount(minCount, Array(length) { count[it] })
        }

        private fun Array<MutableMap<String, Int>>.increment(index: Int, aminoAcid: String) {
            this[index].compute(aminoAcid) { _, u -> (u ?: 0) + 1 }
        }
    }

    operator fun get(locus: Int): Map<String, Int> {
        return count[locus]
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
                lineBuilder.add(base).add(count.toString())
            }

            file.appendText(lineBuilder.toString() + "\n")
        }
    }

}