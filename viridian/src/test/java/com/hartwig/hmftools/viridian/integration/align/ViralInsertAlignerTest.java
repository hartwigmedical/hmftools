package com.hartwig.hmftools.viridian.integration.align;

import static com.hartwig.hmftools.common.bam.CigarUtils.cigarFromStr;
import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.viridian.ViridianTestUtils.noHit;
import static com.hartwig.hmftools.viridian.ViridianTestUtils.samFlags;
import static com.hartwig.hmftools.viridian.ViridianTestUtils.virusReference;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.bwa.MockBwaMemAligner;
import com.hartwig.hmftools.viridian.reference.Virus;
import com.hartwig.hmftools.viridian.reference.VirusGenome;
import com.hartwig.hmftools.viridian.reference.VirusGroup;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.junit.Test;

import htsjdk.samtools.SAMFlag;

public class ViralInsertAlignerTest
{
    private static final VirusGenome GENOME_16 =
            new VirusGenome("hpv16", 7906, new Virus("Human papillomavirus 16", new VirusGroup("HPV 16")));
    private static final VirusGenome GENOME_18 =
            new VirusGenome("hpv18", 7857, new Virus("Human papillomavirus 18", new VirusGroup("HPV 18")));

    private static final int GENOME_16_REF_ID = 0;
    private static final int GENOME_18_REF_ID = 1;

    private static final String SEQUENCE = "ACGTACGTACGTACGTACGTACGTACGTACGT";
    private static final String OTHER_SEQUENCE = "TTGCATTGCATTGCATTGCATTGCATTGCATT";

    @Test
    public void testAlignAllTakesHighestScoringAlignment()
    {
        // The same sequence hits both genomes; only the better-scoring hit is kept.
        List<BwaMemAlignment> hits = List.of(
                forwardAlignment(GENOME_16_REF_ID, 999, 60, 1),
                forwardAlignment(GENOME_18_REF_ID, 499, 45, 4));

        assertEquals(
                new ViralInsertAlignment(GENOME_16, 1000, FORWARD, cigarFromStr("32M"), 60, 1, 32),
                align(SEQUENCE, hits));
    }

    // Near-identical genomes give a sequence equal scores, so the choice must not depend on the order BWA returns them.
    @Test
    public void testAlignAllBreaksScoreTiesDeterministically()
    {
        List<BwaMemAlignment> ascending = List.of(
                forwardAlignment(GENOME_16_REF_ID, 999, 60, 1),
                forwardAlignment(GENOME_18_REF_ID, 499, 60, 1));
        List<BwaMemAlignment> descending = List.of(
                forwardAlignment(GENOME_18_REF_ID, 499, 60, 1),
                forwardAlignment(GENOME_16_REF_ID, 999, 60, 1));

        assertEquals(align(SEQUENCE, ascending), align(SEQUENCE, descending));
        assertEquals(GENOME_16, align(SEQUENCE, descending).genome());
    }

    @Test
    public void testAlignAllReportsReverseStrandAlignment()
    {
        BwaMemAlignment reverse = alignment(samFlags(SAMFlag.READ_REVERSE_STRAND), GENOME_18_REF_ID, 299, 55, 2);

        assertEquals(
                new ViralInsertAlignment(GENOME_18, 300, REVERSE, cigarFromStr("32M"), 55, 2, 32),
                align(SEQUENCE, List.of(reverse)));
    }

    // A sequence matching no genome yields no alignment, and the results stay aligned with the input order.
    @Test
    public void testAlignAllReturnsNullWhereNothingAligned()
    {
        ViralInsertAligner aligner = aligner(Map.of(
                OTHER_SEQUENCE, List.of(noHit()),
                SEQUENCE, List.of(forwardAlignment(GENOME_16_REF_ID, 999, 60, 1))));

        List<ViralInsertAlignment> alignments = aligner.align(List.of(OTHER_SEQUENCE, SEQUENCE));

        assertNull(alignments.get(0));
        assertEquals(GENOME_16, alignments.get(1).genome());
    }

    private static ViralInsertAlignment align(String sequence, List<BwaMemAlignment> hits)
    {
        return aligner(Map.of(sequence, hits)).align(List.of(sequence)).get(0);
    }

    private static ViralInsertAligner aligner(Map<String, List<BwaMemAlignment>> hits)
    {
        return new ViralInsertAligner(MockBwaMemAligner.forQueries(hits), virusReference(GENOME_16, GENOME_18));
    }

    private static BwaMemAlignment forwardAlignment(int refId, int refStart, int score, int editDistance)
    {
        return alignment(0, refId, refStart, score, editDistance);
    }

    // A whole-sequence 32M alignment. BWA reference positions are 0-based.
    private static BwaMemAlignment alignment(int samFlag, int refId, int refStart, int score, int editDistance)
    {
        return new BwaMemAlignment(
                samFlag, refId, refStart, refStart + 32, 0, 32, 60, editDistance, score, 0, "32M", null, null, -1, -1, 0);
    }
}
