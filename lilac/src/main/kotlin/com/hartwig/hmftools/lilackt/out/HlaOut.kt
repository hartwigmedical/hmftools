package com.hartwig.hmftools.lilackt.out

import com.hartwig.hmftools.lilackt.cna.HlaCopyNumber
import com.hartwig.hmftools.lilackt.coverage.HlaAlleleCoverage
import com.hartwig.hmftools.lilackt.coverage.HlaComplexCoverage
import com.hartwig.hmftools.lilackt.variant.SomaticCodingCount
import java.io.File
import java.util.*
import kotlin.math.round

data class HlaOut(val referenceCoverage: HlaComplexCoverage, val tumorCoverage: HlaComplexCoverage, val tumorCopyNumber: List<HlaCopyNumber>, val somaticCodingCount: List<SomaticCodingCount>) {

    companion object {
        fun create(referenceCoverage: HlaComplexCoverage, tumorCoverage: HlaComplexCoverage, tumorCopyNumber: List<HlaCopyNumber>, somaticCodingCount: List<SomaticCodingCount>): HlaOut {
            val sortedCopyNumber = tumorCopyNumber.sortedBy { it.allele }
            val sortedCodingCount = somaticCodingCount.sortedBy { it.allele }

            require(sortedCopyNumber.size == somaticCodingCount.size)
            require(referenceCoverage.alleleCoverage.size == sortedCopyNumber.size)

            return HlaOut(referenceCoverage, tumorCoverage, sortedCopyNumber, sortedCodingCount)
        }
    }

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
                .add("RefTotal")
                .add("RefUnique")
                .add("RefShared")
                .add("RefWild")
                .add("TumorTotal")
                .add("TumorUnique")
                .add("TumorShared")
                .add("TumorWild")
                .add("TumorCopyNumber")
                .add("SomaticMissense")
                .add("SomaticNonsenseOrFrameshift")
                .add("SomaticSplice")
                .add("SomaticSynonymous")
                .add("SomaticInframeIndel")


        return header.toString()
    }

    private fun generateAlleleBody(index: Int): String {
        fun Double.format(digits: Int) = "%.${digits}f".format(this)

        val ref = referenceCoverage.alleleCoverage[index]
        val tumor = if (tumorCoverage.alleleCoverage.isNotEmpty()) tumorCoverage.alleleCoverage[index] else HlaAlleleCoverage(ref.allele, 0, 0.0, 0.0)
        val copyNumber = tumorCopyNumber[index].copyNumber
        val codingCount = somaticCodingCount[index]

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
                .add(copyNumber.format(2))
                .add(codingCount.missense.format(2))
                .add(codingCount.nonsense.format(2))
                .add(codingCount.splice.format(2))
                .add(codingCount.synonymous.format(2))
                .add(codingCount.inframeIndel.format(2))

        return header.toString()
    }


}