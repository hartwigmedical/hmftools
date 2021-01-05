package com.hartwig.hmftools.lilac.nuc

import java.io.File

class NucleotideCount(private val length: Int) {
    val gCount = Array(length) { 0 }
    val aCount = Array(length) { 0 }
    val tCount = Array(length) { 0 }
    val cCount = Array(length) { 0 }

    fun basesAt(index: Int): Collection<Char> {
        val result = mutableSetOf<Char>()
        val depth = depth(index)
        if (gCount[index] > depth / 6) result.add('G')
        if (aCount[index] > depth / 6) result.add('A')
        if (tCount[index] > depth / 6) result.add('T')
        if (cCount[index] > depth / 6) result.add('C')
        return result
    }


    fun depth(index: Int): Int {
        return gCount[index] + aCount[index] + tCount[index] + cCount[index]
    }

    fun isHom(index: Int): Boolean {
        val depth = depth(index)

        val isG = gCount[index] > depth / 6
        val isA = aCount[index] > depth / 6
        val isT = tCount[index] > depth / 6
        val isC = cCount[index] > depth / 6

        var count = 0
        if (isG) count++
        if (isA) count++
        if (isT) count++
        if (isC) count++

        if (count > 1) {
            return false
        }

        return true
    }

    fun writeHorizontally(fileName: String) {
        val file = File(fileName)
        file.writeText("Base\t" + (0 until length).joinToString("\t") + "\n")
        file.appendText("G   \t" + gCount.joinToString("\t") + "\n")
        file.appendText("A   \t" + aCount.joinToString("\t") + "\n")
        file.appendText("T   \t" + tCount.joinToString("\t") + "\n")
        file.appendText("C   \t" + cCount.joinToString("\t") + "\n")
        file.appendText("Hom?\t" + (0 until length).joinToString("\t") { if (isHom(it)) "T" else "F" } + "\n")
    }

    fun writeVertically(fileName: String) {
        val file = File(fileName)
        file.writeText("i\t\tG\t\tA\t\tT\t\tC\t\tHom\n")
        for (i in 0 until length) {
            file.appendText("${i}\t\t${gCount[i]}\t\t${aCount[i]}\t\t${tCount[i]}\t\t${cCount[i]}\t\t${isHom(i)}\n")
        }
    }

}