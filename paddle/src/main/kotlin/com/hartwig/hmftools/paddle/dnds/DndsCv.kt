package com.hartwig.hmftools.paddle.dnds

import com.hartwig.hmftools.paddle.gene.Gene
import java.io.File
import java.nio.file.Files

data class DndsCv(
        val gene: Gene, val nSynonymous: Int,
        val nMissense: Int, val nNonsense: Int, val nSplice: Int, val nIndel: Int,
        val wMissenseCv: Double, val wNonsenseCv: Double, val wSpliceCv: Double, val wIndelCv: Double) {

    companion object {
        fun fromFile(file: String): List<DndsCv> {
            return Files.readAllLines(File(file).toPath()).drop(1).map { fromString(it) }
        }

        private fun fromString(line: String): DndsCv {
            val lineArray = line.split("\t")

            return DndsCv(lineArray[0],
                    lineArray[1].toInt(),
                    lineArray[2].toInt(),
                    lineArray[3].toInt(),
                    lineArray[4].toInt(),
                    lineArray[5].toInt(),
                    lineArray[6].toDouble(),
                    lineArray[7].toDouble(),
                    lineArray[8].toDouble(),
                    lineArray[9].toDouble())
        }

        private fun probability(n: Int, wCv: Double): Double {
            return if (n > 0) {
                ((wCv - 1) / wCv).coerceAtLeast(0.0)
            } else {
                0.0
            }
        }

        private fun expectedDrivers(variantCount: Double, n: Int, wCv: Double): Double {
            return variantCount * probability(n, wCv)
        }
    }

    fun expectedDrivers(impactCount: ImpactCount): ImpactCount {
        if (impactCount.gene != gene) {
            throw IllegalArgumentException("Unexpected gene ${impactCount.gene}")
        }

        return ImpactCount(
                gene,
                expectedMissenseDrivers(impactCount.nMissense),
                expectedNonsenseDrivers(impactCount.nNonsense),
                expectedSpliceDrivers(impactCount.nSplice),
                expectedIndelDrivers(impactCount.nIndel))
    }

    fun expectedMissenseDrivers(nMissense: Double) = expectedDrivers(nMissense, this.nMissense, wMissenseCv)
    fun expectedNonsenseDrivers(nNonsense: Double) = expectedDrivers(nNonsense, this.nNonsense, wNonsenseCv)
    fun expectedSpliceDrivers(nSplice: Double) = expectedDrivers(nSplice, this.nSplice, wSpliceCv)
    fun expectedIndelDrivers(nIndel: Double) = expectedDrivers(nIndel, this.nIndel, wIndelCv)

}