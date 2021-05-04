package com.hartwig.hmftools.lilackt.seq

import com.hartwig.hmftools.lilackt.hla.HlaAllele

data class HlaSequence(val allele: HlaAllele, val rawSequence: String) {
    constructor(allele: String, rawSequence: String) : this(HlaAllele(allele), rawSequence)

    fun copyWithAdditionalSequence(additionalSequence: String): HlaSequence {
        return HlaSequence(allele, rawSequence + additionalSequence)
    }
}