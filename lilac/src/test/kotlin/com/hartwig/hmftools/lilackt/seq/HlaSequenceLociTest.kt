package com.hartwig.hmftools.lilackt.seq

import com.hartwig.hmftools.lilackt.hla.HlaAllele
import org.junit.Assert.*
import org.junit.Test

class HlaSequenceLociTest {

    private val allele = HlaAllele("A*01:01:01:01")
    private val reference = "APRTS..MNV..EPRF"

    @Test
    fun testReference() {
        val victim = HlaSequenceLoci.create(allele, reference, reference)
        assertEquals(12,  victim.sequences.size)
        assertEquals(reference.replace(".", ""), victim.sequence())
        for (sequence in victim.sequences) {
            assertEquals(1, sequence.length)
        }

        assertEquals("APRTSMNVEPRF", victim.sequence())
        assertFalse(victim.containsIndels())
        assertFalse(victim.containsInserts())
        assertFalse(victim.containsDeletes())
    }

    @Test
    fun testInsert() {
        val victimSequence = "-----ET---..----"
        val victim = HlaSequenceLoci.create(allele, victimSequence, reference)
        assertEquals(12,  victim.sequences.size)
        assertEquals("SET",  victim.sequence(4))
        assertEquals("APRTSETMNVEPRF", victim.sequence(0, 11))

        assertTrue(victim.containsIndels())
        assertTrue(victim.containsInserts())
        assertFalse(victim.containsDeletes())

    }

    @Test
    fun testDelete() {
        val victimSequence = ".----...--..---."
        val victim = HlaSequenceLoci.create(allele, victimSequence, reference)
        assertEquals(12,  victim.sequences.size)
        assertEquals(".",  victim.sequence(0))
        assertEquals(".",  victim.sequence(5))
        assertEquals(".",  victim.sequence(11))
        assertEquals(".PRTS.NVEPR.", victim.sequence(0, 11))

        assertTrue(victim.containsIndels())
        assertFalse(victim.containsInserts())
        assertTrue(victim.containsDeletes())
    }

}