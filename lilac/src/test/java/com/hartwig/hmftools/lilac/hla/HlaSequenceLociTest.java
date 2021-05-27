package com.hartwig.hmftools.lilac.hla;

import static com.hartwig.hmftools.lilac.misc.LilacTestUtils.buildTargetIndices;
import static com.hartwig.hmftools.lilac.misc.LilacTestUtils.buildTargetSequence;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertFalse;
import static junit.framework.TestCase.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;
import com.hartwig.hmftools.lilac.seq.SequenceMatchType;

import org.junit.Test;

public class HlaSequenceLociTest
{
    private HlaAllele allele = HlaAllele.fromString("A*01:01:01:01");
    private String reference = "APRTS..MNV..EPRF";

    @Test
    public void testReference()
    {
        HlaSequenceLoci seqLoci = HlaSequenceLoci.create(allele, reference, reference);
        assertEquals(12, seqLoci.getSequences().size());
        assertEquals(reference.replace(".", ""), seqLoci.sequence());

        for(String sequence : seqLoci.getSequences())
        {
            assertEquals(1, sequence.length());
        }

        assertEquals("APRTSMNVEPRF", seqLoci.sequence());
        assertFalse(seqLoci.containsIndels());
        assertFalse(seqLoci.containsInserts());
        assertFalse(seqLoci.containsDeletes());
    }

    @Test
    public void testInsert()
    {
        String sequence = "-----ET---..----";
        HlaSequenceLoci seqLoci = HlaSequenceLoci.create(allele, sequence, reference);
        assertEquals(12, seqLoci.getSequences().size());
        assertEquals("SET", seqLoci.sequence(4));
        assertEquals("APRTSETMNVEPRF", seqLoci.sequence(0, 11));

        assertTrue(seqLoci.containsIndels());
        assertTrue(seqLoci.containsInserts());
        assertFalse(seqLoci.containsDeletes());
    }

    @Test
    public void testDelete()
    {
        String sequence = ".----...--..---.";
        HlaSequenceLoci seqLoci = HlaSequenceLoci.create(allele, sequence, reference);
        assertEquals(12, seqLoci.getSequences().size());
        assertEquals(".", seqLoci.sequence(0));
        assertEquals(".", seqLoci.sequence(5));
        assertEquals(".", seqLoci.sequence(11));
        assertEquals(".PRTS.NVEPR.", seqLoci.sequence(0, 11));

        assertTrue(seqLoci.containsIndels());
        assertFalse(seqLoci.containsInserts());
        assertTrue(seqLoci.containsDeletes());
    }

    @Test
    public void testMatching()
    {
        String refSequence = "ABCDEFGHIJKLMNOP";
        String alleleSequence = "ABXYEFWZIJKLABOP";
        HlaSequenceLoci seqLoci = HlaSequenceLoci.create(allele, alleleSequence, refSequence);

        List<Integer> fragmentLoci = buildTargetIndices(alleleSequence, refSequence);
        String targetSeq = buildTargetSequence(refSequence, fragmentLoci);
        assertEquals(SequenceMatchType.NONE, seqLoci.determineMatchType(targetSeq, fragmentLoci));

        // test against itself
        targetSeq = buildTargetSequence(alleleSequence, fragmentLoci);
        assertEquals(SequenceMatchType.FULL, seqLoci.determineMatchType(targetSeq, fragmentLoci));

        // test an allele with wilcards
        String wildSequence = "ABXYEFW*********";
        HlaSequenceLoci seqLociWild = HlaSequenceLoci.create(allele, wildSequence, refSequence);
        // fragmentLoci = buildTargetIndices(wildSequence, refSequence);

        String fragmentSeq = "ABXYEFWZIJKLABOP";
        fragmentLoci = Lists.newArrayList(5, 6,7,8);
        targetSeq = buildTargetSequence(fragmentSeq, fragmentLoci);
        assertEquals(SequenceMatchType.WILD, seqLociWild.determineMatchType(targetSeq, fragmentLoci));

        // all in the wildcard region
        fragmentLoci = Lists.newArrayList(8,9,10,11);
        targetSeq = buildTargetSequence(fragmentSeq, fragmentLoci);
        assertEquals(SequenceMatchType.WILD, seqLociWild.determineMatchType(targetSeq, fragmentLoci));
    }

}
