package com.hartwig.hmftools.lilac.out

import com.hartwig.hmftools.lilac.coverage.HlaAlleleCoverage
import java.io.File
import java.util.*
import kotlin.math.round

data class HlaOut(val referenceCoverage: List<HlaAlleleCoverage>, val tumorCoverage: List<HlaAlleleCoverage>, val tumorCopyNumber: List<Double>) {

    fun write(fileName: String) {
        val file = File(fileName)
        file.writeText(generateAlleleHeader() + "\n")

        for (i in 0..5) {
            file.appendText(generateAlleleBody(i) + "\n")
        }
    }

    private fun generateAlleleHeader(): String {
        val header = StringJoiner("\t")
                .add("Allele")
                .add("RefTotalCoverage")
                .add("RefUniqueCoverage")
                .add("RefSharedCoverage")
                .add("RefWildCoverage")
                .add("TumorTotalCoverage")
                .add("TumorUniqueCoverage")
                .add("TumorSharedCoverage")
                .add("TumorWildCoverage")
                .add("TumorCopyNumber")

        return header.toString()
    }

    private fun generateAlleleBody(index: Int): String {
        val ref = referenceCoverage[index]
        val tumor = if (tumorCoverage.isNotEmpty()) tumorCoverage[index] else HlaAlleleCoverage(ref.allele, 0, 0.0, 0.0)
        val copyNumber = if (tumorCopyNumber.isNotEmpty()) tumorCopyNumber[index] else 0.0

        val header = StringJoiner("\t")
                .add(ref.allele.asFourDigit().toString())
                .add(round(ref.totalCoverage).toInt().toString())
                .add(ref.uniqueCoverage.toString())
                .add(round(ref.sharedCoverage).toInt().toString())
                .add(round(ref.wildCoverage).toInt().toString())
                .add(round(tumor.totalCoverage).toInt().toString())
                .add(tumor.uniqueCoverage.toString())
                .add(round(tumor.sharedCoverage).toInt().toString())
                .add(round(tumor.wildCoverage).toInt().toString())
                .add(copyNumber.toString())

        return header.toString()
    }


}