package com.hartwig.hmftools.paddle.cohort

import java.io.File
import java.nio.file.Files

data class TumorMutationalLoad(val indel: Int, val mnv: Int, val snv: Int) {
    operator fun plus(other: TumorMutationalLoad): TumorMutationalLoad {
        return TumorMutationalLoad(indel + other.indel, mnv + other.mnv, snv + other.snv)
    }
}

data class TumorMutationalLoadSample(val sample: String, val biallelicLoad: TumorMutationalLoad, val nonBiallelicLoad: TumorMutationalLoad, val totalLoad: TumorMutationalLoad) {

    companion object {
        fun fromFile(file: String): List<TumorMutationalLoadSample> {
            return Files.readAllLines(File(file).toPath()).drop(1).map { fromString(it) }
        }

        private fun fromString(line: String): TumorMutationalLoadSample {
            val lineArray = line.split("\t")

            val biallelicLoad = TumorMutationalLoad(lineArray[1].toInt(), lineArray[2].toInt(), lineArray[3].toInt())
            val nonBiallelicLoad = TumorMutationalLoad(lineArray[4].toInt(), lineArray[5].toInt(), lineArray[6].toInt())
            return TumorMutationalLoadSample(lineArray[0], biallelicLoad, nonBiallelicLoad, biallelicLoad + nonBiallelicLoad)
        }
    }

}
