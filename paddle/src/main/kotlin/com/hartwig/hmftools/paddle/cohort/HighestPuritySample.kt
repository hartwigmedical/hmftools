package com.hartwig.hmftools.paddle.cohort

import com.google.common.collect.Lists
import com.hartwig.hmftools.common.amber.AmberPatient
import com.hartwig.hmftools.common.purple.purity.SamplePurity
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess
import java.io.File
import java.nio.file.Files

data class HighestPuritySample(val sampleId: String, val purity: Double) {

    override fun toString(): String {
        return "$sampleId\t$purity"
    }

    fun toBatchJson(): String {
        return "{\"biopsy\":\"$sampleId\"},"
    }

    companion object {

        fun highestPurityCohort(minPurity: Double, dbAccess: DatabaseAccess): Collection<HighestPuritySample> {
            val purities = dbAccess.readSamplePurityPassingQC(minPurity)
            val patients = dbAccess.readAmberPatients()

            return highestPurityCohort(purities, patients)
        }

        fun highestPurityCohort(purityList: List<SamplePurity>, patientList: List<AmberPatient>): List<HighestPuritySample> {
            val purities = purityList.sortedBy { x -> x.purity() }
            val patients = patientList.associateBy { x -> x.sample() }

            val result = mutableMapOf<Int, HighestPuritySample>()
            for (purity in purities) {
                val sampleId = purity.sampleId()
                if (patients.containsKey(sampleId)) {
                    val patientId = patients[sampleId]!!.patientId()
                    result.put(patientId, HighestPuritySample(sampleId, purity.purity()))
                }
            }

            return result.values.sortedBy { x -> x.sampleId }
        }

        fun readFile(file: String): List<HighestPuritySample> {
            return Files.readAllLines(File(file).toPath()).drop(1).map { fromString(it) }
        }

        private fun fromString(line: String): HighestPuritySample {
            val lineArray = line.split("\t")
            return HighestPuritySample(lineArray[0], lineArray[1].toDouble())
        }

        fun writeFile(filename: String, patients: Collection<HighestPuritySample>) {
            Files.write(File(filename).toPath(), toLines(patients))
        }

        private fun toLines(patients: Collection<HighestPuritySample>): List<String> {
            val lines: MutableList<String> = Lists.newArrayList()
            lines.add("sampleId\tpurity")
            patients.map { it.toString() }.forEach { e: String -> lines.add(e) }
            return lines
        }
    }
}


