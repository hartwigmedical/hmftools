package com.hartwig.hmftools.lilac.hla

import junit.framework.Assert.assertEquals
import org.junit.Test

class HlaReferenceSequenceTest {

    @Test
    fun testRollingKmers() {
        val sequence = HlaReferenceSequence(HlaAllele("A", "01", "01", ""), 6, "NAME", "012345")
        val kmers = sequence.rollingKmers(3)
        assertEquals(4, kmers.size)
    }

}