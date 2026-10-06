package com.hartwig.hmftools.viridian.detection.extract;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.hartwig.hmftools.common.codon.Nucleotides;

import org.junit.Test;

public class ViralKmerIndexTest
{
    private static final String REFERENCE = "ACGATCGGATTCCAAGGTTCAACGTGACCTTAGGCATCGGATACCGTTAC";
    private static final String REFERENCE_WINDOW = REFERENCE.substring(10, 40);

    private static ViralKmerIndex index(String... sequences)
    {
        return ViralKmerIndex.fromSequences(List.of(sequences).stream().map(String::getBytes).toList());
    }

    @Test
    public void testReadSharingReferenceSubsequenceIsDetected()
    {
        assertTrue(index(REFERENCE).hasViralKmer(REFERENCE_WINDOW.getBytes()));
    }

    @Test
    public void testUnrelatedReadIsRejected()
    {
        assertFalse(index(REFERENCE).hasViralKmer("TTTTGGGGCCCCAAAATTTTGGGG".getBytes()));
    }

    @Test
    public void testLowercaseReadMatchesUppercaseReference()
    {
        assertTrue(index(REFERENCE).hasViralKmer(REFERENCE_WINDOW.toLowerCase().getBytes()));
    }

    @Test
    public void testReverseComplementReadIsDetected()
    {
        String reverseComplement = Nucleotides.reverseComplementBases(REFERENCE_WINDOW);
        assertTrue(index(REFERENCE).hasViralKmer(reverseComplement.getBytes()));
    }

    @Test
    public void testReadShorterThanKmerLengthIsRejected()
    {
        assertFalse(index(REFERENCE).hasViralKmer("ACGATCGGATTCCAAG".getBytes()));
    }

    @Test
    public void testKmerSpanningAmbiguousReferenceBaseIsNotIndexed()
    {
        ViralKmerIndex index = index("ACGATCGGANTTCCAAGGT");
        assertEquals(0, index.kmerCount());
        assertFalse(index.hasViralKmer("ACGATCGGATTTCCAAGGT".getBytes()));
    }
}
