package com.hartwig.hmftools.lilackt.hla

import com.hartwig.hmftools.lilackt.nuc.ExpectedAlleles

data class HlaContext(val gene: String, val aminoAcidBoundaries: Set<Int>, val expectedAlleles: ExpectedAlleles)