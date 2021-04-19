package com.hartwig.hmftools.paddle.cohort

import java.io.File
import java.nio.file.Files

data class CohortSample(val sampleId: String, val purity: Double) {

    override fun toString(): String {
        return "$sampleId\t$purity"
    }

    companion object {

        fun readFile(file: String): List<CohortSample> {
            return Files.readAllLines(File(file).toPath()).drop(1).map { fromString(it) }
        }

        private fun fromString(line: String): CohortSample {
            val lineArray = line.split("\t")
            return CohortSample(lineArray[0], lineArray[1].toDouble())
        }
    }
}


