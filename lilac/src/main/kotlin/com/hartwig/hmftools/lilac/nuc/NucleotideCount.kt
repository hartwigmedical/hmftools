package com.hartwig.hmftools.lilac.nuc

import java.io.File

class NucleotideCount(private val length: Int) {
    val gCount = Array(length) { 0 }
    val aCount = Array(length) { 0 }
    val tCount = Array(length) { 0 }
    val cCount = Array(length) { 0 }

    fun depth(index: Int): Int {
        return gCount[index] + aCount[index] + tCount[index] + cCount[index]
    }

    fun isHom(index: Int): Boolean {
        val depth = depth(index)

        val isG = gCount[index] > depth / 4
        val isA = aCount[index] > depth / 4
        val isT = tCount[index] > depth / 4
        val isC = cCount[index] > depth / 4

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

    fun write(fileName: String) {
        val file = File(fileName)
        file.writeText("Base\t" + (0 until length).joinToString("\t") + "\n")
        file.appendText("G   \t" + gCount.joinToString("\t") + "\n")
        file.appendText("A   \t" + aCount.joinToString("\t") + "\n")
        file.appendText("T   \t" + tCount.joinToString("\t") + "\n")
        file.appendText("C   \t" + cCount.joinToString("\t") + "\n")
        file.appendText("Hom?\t" + (0 until length).joinToString("\t") { if (isHom(it)) "T" else "F" } + "\n")

    }

}