package com.hartwig.hmftools.lilac.ext

import htsjdk.samtools.CigarOperator
import htsjdk.samtools.SAMRecord

fun SAMRecord.containsIndel(): Boolean {
    for (cigarElement in this.cigar) {
        if (cigarElement.operator == CigarOperator.I || cigarElement.operator == CigarOperator.D) {
            return true
        }
    }
    return false
}
