package com.hartwig.hmftools.paddle.dnds

import java.io.File
import java.nio.file.Files

data class DndsCv(
        val gene: String,
        val nSynonymous: Int,
        val nMissense: Int,
        val nNonsense: Int,
        val nSplice: Int,
        val nIndel: Int,
        val wMissenseCv: Double,
        val wNonsenseCv: Double,
        val wSpliceCv: Double,
        val wIndelCv: Double) {

    companion object {
        fun fromFile(file: String): List<DndsCv> {
            return Files.readAllLines(File(file).toPath()).drop(1).map { fromString(it) }
        }

        fun fromString(line: String): DndsCv {
            val lineArray = line.split("\t")

            try {
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
            } catch (e: Exception) {
                println(e)
                throw e
            }
        }
    }

}