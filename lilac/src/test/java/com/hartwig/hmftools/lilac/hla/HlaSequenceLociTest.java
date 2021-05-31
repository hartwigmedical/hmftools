package com.hartwig.hmftools.lilac.hla;

import static com.hartwig.hmftools.lilac.misc.LilacTestUtils.buildTargetLoci;
import static com.hartwig.hmftools.lilac.misc.LilacTestUtils.buildTargetSequences;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertFalse;
import static junit.framework.TestCase.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.seq.HlaSequenceFile;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;
import com.hartwig.hmftools.lilac.seq.SequenceMatchType;

import org.junit.Test;

public class HlaSequenceLociTest
{
    private HlaAllele mAllele = HlaAllele.fromString("A*01:01:01:01");
    private HlaAllele mProteinAllele = HlaAllele.fromString("A*01:01");
    private String mReference = "APRTSMNVEPRF";

    @Test
    public void testReference()
    {
        HlaSequenceLoci seqLoci = HlaSequenceLoci.create(mAllele, mReference, mReference);
        assertEquals(12, seqLoci.getSequences().size());
        assertEquals(mReference.replace(".", ""), seqLoci.sequence());

        for(String sequence : seqLoci.getSequences())
        {
            assertEquals(1, sequence.length());
        }

        assertEquals("APRTSMNVEPRF", seqLoci.sequence());
        assertFalse(seqLoci.hasIndels());
        assertFalse(seqLoci.hasInserts());
        assertFalse(seqLoci.hasDeletes());
    }

    @Test
    public void testInsert()
    {
        String sequence = "APRT|SET|MNVEPRF";
        HlaSequenceLoci seqLoci = HlaSequenceFile.createFromReference(mProteinAllele, sequence, true);
        assertEquals(12, seqLoci.getSequences().size());
        assertEquals("SET", seqLoci.sequence(4));
        assertEquals("APRTSETMNVEPRF", seqLoci.sequence(0, 11));

        assertTrue(seqLoci.hasIndels());
        assertTrue(seqLoci.hasInserts());
        assertFalse(seqLoci.hasDeletes());
    }

    @Test
    public void testDelete()
    {
        String sequence = ".PRTS.NVEPR.";
        HlaSequenceLoci seqLoci = HlaSequenceFile.createFromReference(mProteinAllele, sequence, true);
        assertEquals(12, seqLoci.getSequences().size());
        assertEquals(".", seqLoci.sequence(0));
        assertEquals(".", seqLoci.sequence(5));
        assertEquals(".", seqLoci.sequence(11));

        assertTrue(seqLoci.hasIndels());
        assertFalse(seqLoci.hasInserts());
        assertTrue(seqLoci.hasDeletes());
    }

    @Test
    public void testWildcards()
    {
        String sequence = "APRT***MNV**EPRF";
        HlaSequenceLoci seqLoci = HlaSequenceLoci.create(mAllele, sequence, mReference);

        assertFalse(seqLoci.hasInserts());
        assertFalse(seqLoci.hasDeletes());
        assertTrue(seqLoci.hasWildcards());

        seqLoci.setExonBoundaryWildcardsWildcards(Lists.newArrayList(2));
        assertFalse(seqLoci.hasExonBoundaryWildcards());

        seqLoci.setExonBoundaryWildcardsWildcards(Lists.newArrayList(5));
        assertTrue(seqLoci.hasExonBoundaryWildcards());
    }

    @Test
    public void testMatching()
    {
        String refSequence = "ABCDEFGHIJKLMNOP";
        String alleleSequence = "ABXYEFWZIJKLABOP";
        HlaSequenceLoci seqLoci = HlaSequenceLoci.create(mAllele, alleleSequence, refSequence);

        List<Integer> fragmentLoci = buildTargetLoci(alleleSequence, refSequence);
        List<String> targetSeq = buildTargetSequences(refSequence, fragmentLoci);
        assertEquals(SequenceMatchType.NONE, seqLoci.determineMatchType(targetSeq, fragmentLoci));

        // test against itself
        targetSeq = buildTargetSequences(alleleSequence, fragmentLoci);
        assertEquals(SequenceMatchType.FULL, seqLoci.determineMatchType(targetSeq, fragmentLoci));

        // test an allele with wilcards
        String wildSequence = "ABXYEFW*********";
        HlaSequenceLoci seqLociWild = HlaSequenceLoci.create(mAllele, wildSequence, refSequence);
        // fragmentLoci = buildTargetIndices(wildSequence, refSequence);

        String fragmentSeq = "ABXYEFWZIJKLABOP";
        fragmentLoci = Lists.newArrayList(5, 6,7,8);
        targetSeq = buildTargetSequences(fragmentSeq, fragmentLoci);
        assertEquals(SequenceMatchType.WILD, seqLociWild.determineMatchType(targetSeq, fragmentLoci));

        // all in the wildcard region
        fragmentLoci = Lists.newArrayList(8,9,10,11);
        targetSeq = buildTargetSequences(fragmentSeq, fragmentLoci);
        assertEquals(SequenceMatchType.WILD, seqLociWild.determineMatchType(targetSeq, fragmentLoci));
    }

}
