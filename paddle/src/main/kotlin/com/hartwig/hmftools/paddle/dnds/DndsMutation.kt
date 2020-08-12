package com.hartwig.hmftools.paddle.dnds

import java.io.File
import java.nio.file.Files

data class DndsMutation(
        val sample: String,
        val contig: String, val position: Int, val ref: String, val alt: String,
        val worstCodingEffect: String, val canonicalCodingEffect: String, val repeatCount: Int, val biallelic: Boolean, val hotspot: Boolean,
        val gene: String, val dndsImpact: String) {

    companion object {
        fun fromFile(file: String): List<DndsMutation> {
            return Files.readAllLines(File(file).toPath()).drop(1).map { fromString(it) }
        }

        fun fromString(line: String): DndsMutation {
            val array = line.split("\t")
            return DndsMutation(
                    array[0],
                    array[1],
                    array[2].toInt(),
                    array[3],
                    array[4],
                    array[5],
                    array[6],
                    array[7].toInt(),
                    array[8].toBoolean(),
                    array[9] == "HOTSPOT",
                    array[10],
                    array[11])
        }


    }


}
