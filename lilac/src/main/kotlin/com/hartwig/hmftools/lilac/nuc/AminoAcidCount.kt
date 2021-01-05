package com.hartwig.hmftools.lilac.nuc

import java.io.File
import java.util.*
import kotlin.math.min

class AminoAcidCount(private val length: Int) {
    private val count = Array(length) { mutableMapOf<Char, Int>() }

    fun increment(index: Int, aminoAcid: Char) {
        if (aminoAcid != '.') {
            count[index].compute(aminoAcid) { _, u -> (u ?: 0) + 1 }
        }
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