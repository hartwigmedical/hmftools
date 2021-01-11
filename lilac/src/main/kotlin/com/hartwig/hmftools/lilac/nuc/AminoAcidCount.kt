package com.hartwig.hmftools.lilac.nuc

import com.hartwig.hmftools.lilac.read.Fragment
import java.io.File
import java.util.*
import kotlin.math.min

class AminoAcidCount(val length: Int) {

    companion object {
        operator fun invoke(minQual: Int, fragments: List<Fragment>): AminoAcidCount {
            val length = fragments.map { it.aminoAcidIndices().max() ?: -1 }.max()!! + 1
            val result = AminoAcidCount(length)

            for (fragment in fragments) {
                for (index in fragment.aminoAcidIndices()) {
                    val aminoAcid = fragment.aminoAcid(index, minQual)
                    if (aminoAcid != '.') {
                        result.increment(index, aminoAcid)
                    }
                }
            }
            return result
        }
    }

    private val count = Array(length) { mutableMapOf<Char, Int>() }

    private fun increment(index: Int, aminoAcid: Char) {
        if (aminoAcid != '.') {
            count[index].compute(aminoAcid) { _, u -> (u ?: 0) + 1 }
        }
    }

    fun heterozygousIndices(minCount: Int): List<Int> {
        val result = mutableListOf<Int>()
        for (i in count.indices) {
            if (isHeterozygous(minCount, i)) {
                result.add(i)
            }
        }

        return result
    }

    private fun isHeterozygous(minCount: Int, index: Int): Boolean {
        val mapAtIndex = count[index].filter { it.value >= minCount }
        return mapAtIndex.size > 1
    }

    fun aminoAcidAt(index: Int, minCount: Int = 2): Collection<Char> {
        val result = mutableSetOf<Char>()
        val indexMap = count[index]
        for ((aa, count) in indexMap) {
            if (aa != '.' && count >= minCount) {
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