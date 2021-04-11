package com.hartwig.hmftools.lilac.qc

import java.io.File


data class LilacQC(val status: LilacQCStatus, val aminoAcidQC: AminoAcidQC, val bamQC: BamQC, val coverageQC: CoverageQC, val haplotypeQC: HaplotypeQC, val somaticVariantQC: SomaticVariantQC) {

    companion object {
        fun create(aminoAcidQC: AminoAcidQC, bamQC: BamQC, coverageQC: CoverageQC, haplotypeQC: HaplotypeQC, somaticVariantQC: SomaticVariantQC): LilacQC {
            var status = LilacQCStatus.PASS

            if (haplotypeQC.unusedHaplotypes > 0) {
                status = LilacQCStatus.WARN_UNMATCHED_HAPLOTYPE
            }

            if (coverageQC.aTypes == 0 || coverageQC.bTypes == 0 || coverageQC.cTypes == 0 ) {
                status = LilacQCStatus.WARN_UNMATCHED_TYPE
            }

            if (bamQC.discardedIndelFragments > 0) {
                status = LilacQCStatus.WARN_UNMATCHED_INDEL
            }

            if (somaticVariantQC.unmatchedVariants()) {
                status = LilacQCStatus.WARN_UNMATCHED_SOMATIC_VARIANT
            }

            if (coverageQC.percentWildcard > 0) {
                status = LilacQCStatus.WARN_WILDCARD_MATCH
            }

            return LilacQC(status, aminoAcidQC, bamQC, coverageQC, haplotypeQC, somaticVariantQC)
        }
    }

    fun header(): List<String> {
        return listOf("status") + aminoAcidQC.header() + bamQC.header() + coverageQC.header() + haplotypeQC.header() + somaticVariantQC.header()
    }

    fun body(): List<String> {
        return listOf(status.toString()) + aminoAcidQC.body() + bamQC.body() + coverageQC.body() + haplotypeQC.body() + somaticVariantQC.body()
    }

    fun writefile(fileName: String) {
        val header = header().joinToString(separator = "\t")
        val body = body().joinToString(separator = "\t")

        val file = File(fileName)
        file.writeText(header + "\n")
        file.appendText(body + "\n")
    }
}

