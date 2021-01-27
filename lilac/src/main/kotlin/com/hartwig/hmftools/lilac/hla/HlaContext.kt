package com.hartwig.hmftools.lilac.hla

import com.hartwig.hmftools.lilac.nuc.ExpectedAlleles

data class HlaContext(val gene: String, val aminoAcidBoundaries: Set<Int>, val expectedAlleles: ExpectedAlleles) {

    companion object {
        fun hlaA(aBoundaries: Set<Int>, bBoundaries: Set<Int>, cBoundaries: Set<Int>): HlaContext {
            return HlaContext("A", aBoundaries, ExpectedAlleles.create(aBoundaries, bBoundaries, cBoundaries))
        }

        fun hlaB(aBoundaries: Set<Int>, bBoundaries: Set<Int>, cBoundaries: Set<Int>): HlaContext {
            return HlaContext("B", bBoundaries, ExpectedAlleles.create(bBoundaries, aBoundaries, cBoundaries))
        }

        fun hlaC(aBoundaries: Set<Int>, bBoundaries: Set<Int>, cBoundaries: Set<Int>): HlaContext {
            return HlaContext("C", cBoundaries, ExpectedAlleles.create(cBoundaries, aBoundaries, bBoundaries))
        }

    }


}