package com.hartwig.hmftools.lilac.hla

import com.hartwig.hmftools.lilac.amino.AminoAcidFragment
import com.hartwig.hmftools.lilac.nuc.ExpectedAlleles

data class HlaContext(val gene: String, val aminoAcidBoundaries: Set<Int>, val expectedAlleles: ExpectedAlleles, val fragments: List<AminoAcidFragment>) {

    companion object {
        fun hlaA(aBoundaries: Set<Int>, bBoundaries: Set<Int>, cBoundaries: Set<Int>, fragments: List<AminoAcidFragment>): HlaContext {
            return HlaContext("A", aBoundaries, ExpectedAlleles.create(aBoundaries, bBoundaries, cBoundaries), fragments)
        }

        fun hlaB(aBoundaries: Set<Int>, bBoundaries: Set<Int>, cBoundaries: Set<Int>, fragments: List<AminoAcidFragment>): HlaContext {
            return HlaContext("B", bBoundaries, ExpectedAlleles.create(bBoundaries, aBoundaries, cBoundaries), fragments)
        }

        fun hlaC(aBoundaries: Set<Int>, bBoundaries: Set<Int>, cBoundaries: Set<Int>, fragments: List<AminoAcidFragment>): HlaContext {
            return HlaContext("C", cBoundaries, ExpectedAlleles.create(cBoundaries, aBoundaries, bBoundaries), fragments)
        }

    }


}