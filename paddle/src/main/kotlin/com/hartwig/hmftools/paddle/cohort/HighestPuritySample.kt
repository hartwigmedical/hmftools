package com.hartwig.hmftools.paddle.cohort

import com.google.common.collect.Lists
import com.hartwig.hmftools.common.amber.AmberPatient
import com.hartwig.hmftools.common.purple.purity.SamplePurity
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess
import java.io.File
import java.nio.file.Files

data class HighestPuritySample(val patientId: Int, val sampleId: String, val purity: Double) {

    override fun toString(): String {
        return "$patientId\t$sampleId\t$purity"
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
                    result.put(patientId, HighestPuritySample(patientId, sampleId, purity.purity()))
                }
            }

            return result.values.sortedBy { x -> x.patientId }
        }

        fun writeFile(filename: String, patients: Collection<HighestPuritySample>) {
            Files.write(File(filename).toPath(), toLines(patients))
        }

        private fun toLines(patients: Collection<HighestPuritySample>): List<String> {
            val lines: MutableList<String> = Lists.newArrayList()
            lines.add("patientId\tsampleId\tpurity")
            patients.map { it.toString() }.forEach { e: String -> lines.add(e) }
            return lines
        }
    }
}


