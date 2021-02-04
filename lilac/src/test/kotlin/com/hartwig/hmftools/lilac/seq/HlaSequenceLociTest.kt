package com.hartwig.hmftools.lilac.seq

import com.hartwig.hmftools.lilac.hla.HlaAllele
import org.junit.Assert.assertEquals
import org.junit.Test

class HlaSequenceLociTest {

    val allele = HlaAllele("A*01:01:01:01")
    val reference = "APRTS..MNV..EPRF"

    @Test
    fun testReference() {
        val victim = HlaSequenceLoci.create(allele, reference, reference)
        assertEquals(12,  victim.sequences.size)
        assertEquals(reference.replace(".", ""), victim.sequence())
        for (sequence in victim.sequences) {
            assertEquals(1, sequence.length)
        }

        assertEquals("APRTSMNVEPRF", victim.sequence())
    }

    @Test
    fun testInsert() {
        val victimSequence = "-----ET---..----"
        val victim = HlaSequenceLoci.create(allele, victimSequence, reference)
        assertEquals(12,  victim.sequences.size)
        assertEquals("SET",  victim.sequence(4))
        assertEquals("APRTSETMNVEPRF", victim.sequence())
    }

    @Test
    fun testDelete() {
        val victimSequence = ".----...--..---."
        val victim = HlaSequenceLoci.create(allele, victimSequence, reference)
        assertEquals(12,  victim.sequences.size)
        assertEquals(".",  victim.sequence(0))
        assertEquals(".",  victim.sequence(5))
        assertEquals(".",  victim.sequence(11))
        assertEquals("PRTSNVEPR", victim.sequence())
    }

}