package com.hartwig.hmftools.lilac.qc

import java.io.File


data class LilacQC(val aminoAcidQC: AminoAcidQC, val bamQC: BamQC, val coverageQC: CoverageQC, val haplotypeQC: HaplotypeQC, val somaticVariantQC: SomaticVariantQC) {

    fun header(): List<String> {
        return aminoAcidQC.header() + bamQC.header() + coverageQC.header() + haplotypeQC.header() + somaticVariantQC.header()
    }

    fun body(): List<String> {
        return aminoAcidQC.body() + bamQC.body() + coverageQC.body() + haplotypeQC.body() + somaticVariantQC.body()
    }

    fun writefile(fileName: String) {
        val header = header().joinToString(separator = "\t")
        val body = body().joinToString(separator = "\t")

        val file = File(fileName)
        file.writeText(header + "\n")
        file.appendText(body + "\n")
    }
}

