package com.hartwig.hmftools.lilackt.hla

import com.hartwig.hmftools.lilackt.nuc.ExpectedAlleles
import com.hartwig.hmftools.lilackt.nuc.NucleotideGeneEnrichment

class HlaContextFactory(private val aBoundaries: Set<Int>, private val bBoundaries: Set<Int>, private val cBoundaries: Set<Int>) {
    private val nucleotideGeneEnrichment = NucleotideGeneEnrichment(aBoundaries, bBoundaries, cBoundaries)

    fun hlaA(): HlaContext {
        val expectedAlleles = ExpectedAlleles.expectedAlleles(nucleotideGeneEnrichment.aFilterB, nucleotideGeneEnrichment.aFilterC)
        return HlaContext("A", aBoundaries, expectedAlleles)
    }

    fun hlaB(): HlaContext {
        val expectedAlleles = ExpectedAlleles.expectedAlleles(nucleotideGeneEnrichment.bFilterA, nucleotideGeneEnrichment.bFilterC)
        return HlaContext("B", bBoundaries, expectedAlleles)
    }

    fun hlaC(): HlaContext {
        val expectedAlleles = ExpectedAlleles.expectedAlleles(nucleotideGeneEnrichment.cFilterA, nucleotideGeneEnrichment.cFilterB)
        return HlaContext("C", cBoundaries, expectedAlleles)
    }

}