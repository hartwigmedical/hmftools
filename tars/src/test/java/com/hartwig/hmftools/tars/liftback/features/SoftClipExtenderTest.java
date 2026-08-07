package com.hartwig.hmftools.tars.liftback.features;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.bases;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.primaryRecord;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.selfAlignment;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.hartwig.hmftools.tars.common.BwaScoring;
import com.hartwig.hmftools.tars.liftback.LiftedAlignment;
import com.hartwig.hmftools.tars.liftback.LiftedRecord;
import com.hartwig.hmftools.tars.liftback.TarsTestFixtures.TestGenome;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;

public class SoftClipExtenderTest
{
    private static TestGenome genome()
    {
        return new TestGenome().with(CHR_1, 5000, 'A');
    }

    private static SoftClipExtender extender(final TestGenome genome)
    {
        return new SoftClipExtender(genome.asRefGenome());
    }

    private static LiftedAlignment refAlt(final String chrom, final int pos, final String cigar, final boolean forwardStrand)
    {
        return new LiftedAlignment(chrom, pos, cigar, 0, false, false, forwardStrand, 0);
    }

    @Test
    public void testExtendsTerminalSoftclipCases()
    {
        SoftClipExtender trailingExtender = extender(genome().set(CHR_1, 120, "CCCCC"));
        SoftClipExtender leadingExtender = extender(genome().set(CHR_1, 105, "CCCCC"));
        SoftClipExtender mismatchExtender = extender(genome());

        LiftedAlignment trailing = trailingExtender.extend(
                selfAlignment(CHR_1, 100, "20M5S"), bases("A".repeat(20) + "C".repeat(5)));
        LiftedAlignment leading = leadingExtender.extend(
                selfAlignment(CHR_1, 110, "5S20M"), bases("C".repeat(5) + "A".repeat(20)));
        LiftedAlignment mismatch = mismatchExtender.extend(
                selfAlignment(CHR_1, 100, "20M5S"), bases("A".repeat(20) + "C".repeat(5)));

        // Trailing soft clip matches the next reference bases.
        assertEquals(100, trailing.LiftedPos);
        assertEquals("25M", trailing.LiftedCigar);

        // Leading soft clip matches the previous reference bases.
        assertEquals(105, leading.LiftedPos);
        assertEquals("25M", leading.LiftedCigar);

        // Mismatching clipped bases are left alone.
        assertEquals(100, mismatch.LiftedPos);
        assertEquals("20M5S", mismatch.LiftedCigar);
    }

    @Test
    public void testAppliesToEverySurvivingCandidate()
    {
        TestGenome genome = genome()
                .set(CHR_1, 120, "CCCCC")
                .set(CHR_1, 295, "CCCCC")
                .set(CHR_1, 420, "CCCCC");
        SoftClipExtender extender = extender(genome);

        LiftedAlignment self = selfAlignment(CHR_1, 100, "20M5S");
        LiftedAlignment alt = refAlt(CHR_1, 300, "5S20M", true);
        LiftedAlignment dropped = refAlt(CHR_1, 400, "20M5S", true);
        dropped.Dropped = true;
        List<LiftedAlignment> alignments = new ArrayList<>(List.of(self, alt, dropped));

        SAMRecord record = primaryRecord(CHR_1, 100, "20M5S");
        record.setReadBases(bases("C".repeat(5) + "A".repeat(15) + "C".repeat(5)));

        LiftedRecord liftedRecord = new LiftedRecord(60, 1, "", 0, alignments)
                .withExtendedSoftClips(extender, record);
        alignments = liftedRecord.liftedAlignments();

        // Every surviving candidate uses the same feature pass.
        assertEquals("25M", alignments.get(0).LiftedCigar);
        assertEquals("25M", alignments.get(1).LiftedCigar);
        assertEquals(295, alignments.get(1).LiftedPos);

        // Dropped candidates stay untouched.
        assertEquals("20M5S", alignments.get(2).LiftedCigar);
        assertTrue(alignments.get(2).Dropped);
    }

    @Test
    public void testScoresCandidatesAfterExtension()
    {
        TestGenome genome = genome().set(CHR_1, 2000, "T".repeat(50));
        SoftClipExtender extender = extender(genome);

        LiftedAlignment self = selfAlignment(CHR_1, 100, "50M");
        LiftedAlignment oppositeStrandAlt = refAlt(CHR_1, 2000, "50M", false);
        List<LiftedAlignment> alignments = new ArrayList<>(List.of(self, oppositeStrandAlt));

        SAMRecord record = primaryRecord(CHR_1, 100, "50M");
        byte[] forward = bases("A".repeat(50));
        record.setReadBases(forward);

        extender.scoreCandidates(alignments, record);

        byte[] reverseComplement = Arrays.copyOf(forward, forward.length);
        SequenceUtil.reverseComplement(reverseComplement);

        // Surviving candidates are scored on their own strand.
        assertEquals(50, self.GenomicScore);
        assertEquals(BwaScoring.genomicScore(genome.asRefGenome(), CHR_1, 2000, "50M", reverseComplement),
                oppositeStrandAlt.GenomicScore);
        assertTrue(BwaScoring.genomicScore(genome.asRefGenome(), CHR_1, 2000, "50M", forward)
                < oppositeStrandAlt.GenomicScore);
    }

    @Test
    public void testSupplementaryTagDoesNotSuppressCandidateScoring()
    {
        SoftClipExtender extender = extender(genome());

        List<LiftedAlignment> alignments = new ArrayList<>(
                List.of(selfAlignment(CHR_1, 100, "50M"), refAlt(CHR_1, 2000, "50M", true)));
        SAMRecord record = primaryRecord(CHR_1, 100, "50M");
        record.setReadBases(bases("A".repeat(50)));
        record.setAttribute("SA", CHR_1 + ",2000,+,50M,0,0;");

        extender.scoreCandidates(alignments, record);

        assertEquals(50, alignments.get(0).GenomicScore);
        assertEquals(50, alignments.get(1).GenomicScore);
    }
}
