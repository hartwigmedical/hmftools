package com.hartwig.hmftools.lilac.hla;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertFalse;
import static junit.framework.TestCase.assertTrue;

import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;

import org.junit.Test;

public class HlaSequenceLociTest
{
    private HlaAllele allele = HlaAllele.fromString("A*01:01:01:01");
    private String reference = "APRTS..MNV..EPRF";

    @Test
    public void testReference()
    {
        HlaSequenceLoci victim = HlaSequenceLoci.create(allele, reference, reference);
        assertEquals(12, victim.getSequences().size());
        assertEquals(reference.replace(".", ""), victim.sequence());

        for(String sequence : victim.getSequences())
        {
            assertEquals(1, sequence.length());
        }

        assertEquals("APRTSMNVEPRF", victim.sequence());
        assertFalse(victim.containsIndels());
        assertFalse(victim.containsInserts());
        assertFalse(victim.containsDeletes());
    }

    @Test
    public void testInsert()
    {
        String victimSequence = "-----ET---..----";
        HlaSequenceLoci victim = HlaSequenceLoci.create(allele, victimSequence, reference);
        assertEquals(12, victim.getSequences().size());
        assertEquals("SET", victim.sequence(4));
        assertEquals("APRTSETMNVEPRF", victim.sequence(0, 11));

        assertTrue(victim.containsIndels());
        assertTrue(victim.containsInserts());
        assertFalse(victim.containsDeletes());
    }

    @Test
    public void testDelete()
    {
        String victimSequence = ".----...--..---.";
        HlaSequenceLoci victim = HlaSequenceLoci.create(allele, victimSequence, reference);
        assertEquals(12, victim.getSequences().size());
        assertEquals(".", victim.sequence(0));
        assertEquals(".", victim.sequence(5));
        assertEquals(".", victim.sequence(11));
        assertEquals(".PRTS.NVEPR.", victim.sequence(0, 11));

        assertTrue(victim.containsIndels());
        assertFalse(victim.containsInserts());
        assertTrue(victim.containsDeletes());
    }

}
