package com.hartwig.hmftools.paddle.dnds

import com.hartwig.hmftools.paddle.Gene
import java.io.File
import java.nio.file.Files

data class DndsCv(val n: Int, val wCv: Double) {
    fun expectedDrivers(variantCount: Int): Double {
        return variantCount * probability(n, wCv)
    }

    private fun probability(n: Int, wCv: Double): Double {
        return if (n > 0) {
            ((wCv - 1) / wCv).coerceAtLeast(0.0)
        } else {
            0.0
        }
    }
}

data class DndsCvGene(val gene: Gene, val nSynonymous: Int,
                      val missense: DndsCv, val nonsense: DndsCv, val splice: DndsCv, val indel: DndsCv) {

    companion object {
        fun fromFile(file: String): List<DndsCvGene> {
            return Files.readAllLines(File(file).toPath()).drop(1).map { fromString(it) }
        }

        private fun fromString(line: String): DndsCvGene {
            val lineArray = line.split("\t")

            return DndsCvGene(lineArray[0], lineArray[1].toInt(),
                    DndsCv(lineArray[2].toInt(), lineArray[6].toDouble()),
                    DndsCv(lineArray[3].toInt(), lineArray[7].toDouble()),
                    DndsCv(lineArray[4].toInt(), lineArray[8].toDouble()),
                    DndsCv(lineArray[5].toInt(), lineArray[9].toDouble()))
        }
    }
}