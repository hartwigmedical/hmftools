package com.hartwig.hmftools.chord;

import static org.junit.Assert.assertEquals;

import com.google.common.io.Resources;
import com.hartwig.hmftools.chord.indel.IndelVariant;
import com.hartwig.hmftools.chord.common.SmallVariant;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;

import org.junit.Test;

public class IndelVariantTest
{
    private static final String REF_GENOME_FILE = Resources.getResource("fasta/indel_regions.fasta").getPath();
    private static final RefGenomeSource REF_GENOME = RefGenomeSource.loadRefGenome(REF_GENOME_FILE);

    @Test
    public void canGetDeletionFlankSequences()
    {
        // Based on deletion at chr:222493574 in COLO829v003T
        IndelVariant deletion = new IndelVariant(new SmallVariant("chr1", 21, "CAA", "C"));

        String leftFlank = deletion.getFlankingRefBasesLeft(REF_GENOME, 1);
        String rightFlank = deletion.getFlankingRefBasesRight(REF_GENOME, 3);

        System.out.println("Left flank: " + leftFlank);
        System.out.println("Right flank: " + rightFlank);

        assertEquals(deletion.mIndelLength, 2);
        assertEquals("CC", leftFlank);
        assertEquals("AAGTCT", rightFlank);
    }

    @Test
    public void canGetInsertionFlankSequences()
    {
        // Based on insertion at chr17:79329838 in COLO829v003T
        IndelVariant insertion = new IndelVariant(new SmallVariant("chr17", 21, "G", "GGGTGT"));

        String leftFlank = insertion.getFlankingRefBasesLeft(REF_GENOME, 1);
        String rightFlank = insertion.getFlankingRefBasesRight(REF_GENOME, 3);

        System.out.println("Left flank: " + leftFlank);
        System.out.println("Right flank: " + rightFlank);

        assertEquals(insertion.mIndelLength, 5);
        assertEquals("GGTCG", leftFlank);
        assertEquals("GGTGTGGTGGCTCAC", rightFlank);
    }

    @Test
    public void canCountHomoglousBases()
    {
        String indelSequence = "AAGGCACTGGAAAGAATATAGAGA";
        String flankSequence = "AAGGCACTCTCCACTCTGCTTTGTACATTTCTCATTACTGTATTAGCATCTTATAATCATAAATAACACCTT";

        int homologousBases = IndelVariant.countHomologousBases(indelSequence, flankSequence);

        assertEquals(8, homologousBases);
    }

    @Test
    public void canCountRepeatUnits()
    {
        String indelSequence = "CT";
        String rightFlankSequence = "CTCTCAAA";

        int repeatUnits = IndelVariant.countRepeatUnits(indelSequence, rightFlankSequence);

        assertEquals(3, repeatUnits);
    }

    @Test
    public void cantGetIndelSequence()
    {
        IndelVariant deletion = new IndelVariant(new SmallVariant("chr1", 1, "AGGG", "A"));
        assertEquals("GGG", deletion.getIndelSequence());

        IndelVariant insertion = new IndelVariant(new SmallVariant("chr1", 1, "C", "CTTT"));
        assertEquals("TTT", insertion.getIndelSequence());
    };
}
