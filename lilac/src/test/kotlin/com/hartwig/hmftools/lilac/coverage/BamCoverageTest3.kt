package com.hartwig.hmftools.lilac.coverage

import com.hartwig.hmftools.common.codon.Codons
import org.junit.Test

class BamCoverageTest3 {

    companion object {
        const val EXON1 = "SHSMRYFDTAVSRPGRGEPRFISVGYVDDTQFVRFDSDAASPRGEPRAPWVEQEGPEYWDRETQNYKRQAQADRVSLRNLRGYYNQSED"
        const val EXON2 = "PPKTHVTHHPLSDHEATLRCWALGFYPAEITLTWQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGQEQRYTCHMQHEGLQEPLTLSW"
    }

    @Test
    fun testProcessDna() {
        val exon1Coverage = ExonCoverage(EXON1)
        val exon2Coverage = ExonCoverage(EXON2)
        val victim = BamCoverage3(3, listOf(exon1Coverage, exon2Coverage))

        victim.processDna(("HSMRY").toDna())
        println(exon1Coverage.coverage.joinToString())
    }

    private fun String.startCodons(codons: Int): String {
        return this.substring(0, codons).toDna()
    }

    private fun String.endCodons(codons: Int): String {
        return this.substring(this.length - codons, this.length).toDna()
    }

    private fun String.toDna(): String {
        return Codons.codon(this)
    }

}