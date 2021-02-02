package com.hartwig.hmftools.lilac.hla

import com.hartwig.hmftools.lilac.nuc.ExpectedAlleles

data class HlaContext(val gene: String, val aminoAcidBoundaries: Set<Int>, val expectedAlleles: ExpectedAlleles)