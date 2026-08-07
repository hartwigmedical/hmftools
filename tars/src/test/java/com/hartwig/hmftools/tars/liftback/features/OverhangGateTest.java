package com.hartwig.hmftools.tars.liftback.features;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.bases;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.primaryRecord;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.selfAlignment;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.tars.liftback.LiftedAlignment;
import com.hartwig.hmftools.tars.liftback.TarsTestFixtures.TestGenome;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

// README Step 1: short terminal overhangs (<= 12M) are scored against the reference. With soft clip, keep only if
// overhang AS > 5. With multiple junctions and no soft clip, keep if AS > 0; otherwise collapse only when read-through
// scores better. A single junction with no soft clip is left alone.
public class OverhangGateTest
{
    private static OverhangGate gate(final TestGenome genome)
    {
        return new OverhangGate(genome.asRefGenome());
    }

    private static TestGenome genome()
    {
        return new TestGenome().with(CHR_1, 5000, 'A');
    }

    private static LiftedAlignment refAlt(final String chrom, final int pos, final String cigar, final boolean forwardStrand)
    {
        return new LiftedAlignment(chrom, pos, cigar, 0, false, false, forwardStrand, 0);
    }

    @Test
    public void testSoftclipAdjacentJunctionCases()
    {
        OverhangGate gate = gate(genome().set(CHR_1, 121, "CCCCCC"));

        OverhangGate.Outcome supported = gate.gate(
                selfAlignment(CHR_1, 1, "20M100N6M4S"), bases("C".repeat(30)));
        OverhangGate.Outcome unsupported = gate(genome())
                .gate(selfAlignment(CHR_1, 1000, "100M83N3M48S"), bases("C".repeat(151)));
        OverhangGate.Outcome longAnchor = gate(genome())
                .gate(selfAlignment(CHR_1, 1, "20M100N15M4S"), bases("C".repeat(39)));
        OverhangGate.Outcome noSoftclip = gate(genome())
                .gate(selfAlignment(CHR_1, 1, "20M100N4M"), bases("C".repeat(24)));
        OverhangGate.Outcome leading = gate(genome().set(CHR_1, 97, "CCCCCCCCC"))
                .gate(selfAlignment(CHR_1, 1, "4S5M100N20M"), bases("A".repeat(29)));

        // 6M scores above the soft-clip threshold, so the junction is kept.
        assertEquals("20M100N6M4S", supported.alignment().LiftedCigar);

        // 3M is below the threshold, so the N folds into the terminal soft clip.
        assertEquals("100M51S", unsupported.alignment().LiftedCigar);
        assertTrue(unsupported.dropped());

        // Long terminal anchors are trusted without scoring.
        assertEquals("20M100N15M4S", longAnchor.alignment().LiftedCigar);

        // A single junction with no adjacent soft clip is bwa's call and is left alone.
        assertEquals("20M100N4M", noSoftclip.alignment().LiftedCigar);

        // Leading collapse moves the start to the near exon and converts the overhang to soft clip.
        assertEquals("9S20M", leading.alignment().LiftedCigar);
        assertEquals(106, leading.alignment().LiftedPos);
        assertTrue(leading.dropped());
    }

    @Test
    public void testMultiJunctionReadThroughCases()
    {
        OverhangGate.Outcome positiveAnchor = gate(genome())
                .gate(selfAlignment(CHR_1, 1, "20M50N5M60N4M"), bases("A".repeat(29)));
        OverhangGate.Outcome readThroughNotBetter = gate(genome())
                .gate(selfAlignment(CHR_1, 1, "20M50N5M60N4M"), bases("C".repeat(29)));
        OverhangGate.Outcome readThroughBetter = gate(genome().set(CHR_1, 136, "CCCC"))
                .gate(selfAlignment(CHR_1, 1, "20M50N5M60N4M"), bases("A".repeat(29)));
        OverhangGate.Outcome terminalKept = gate(genome())
                .gate(selfAlignment(CHR_1, 1, "20M50N4M60N40M"), bases("A".repeat(64)));

        // Positive terminal-anchor score keeps a no-softclip multi-junction alignment.
        assertEquals("20M50N5M60N4M", positiveAnchor.alignment().LiftedCigar);

        // If read-through is no better than the anchor, the junction is kept.
        assertEquals("20M50N5M60N4M", readThroughNotBetter.alignment().LiftedCigar);

        // If read-through scores better, only the terminal N collapses; the remaining N means the placement is kept.
        assertEquals("20M50N9M", readThroughBetter.alignment().LiftedCigar);
        assertFalse(readThroughBetter.dropped());

        // Once the terminal junction is kept, interior junctions behind it are not scored.
        assertEquals("20M50N4M60N40M", terminalKept.alignment().LiftedCigar);
    }

    @Test
    public void testAppliesGateToCandidateList()
    {
        OverhangGate gate = gate(genome());

        List<LiftedAlignment> selfOnly = new ArrayList<>(List.of(selfAlignment(CHR_1, 1, "20M100N8M4S")));
        SAMRecord selfRecord = primaryRecord(CHR_1, 1, "20M100N8M4S");
        selfRecord.setReadBases(bases("C".repeat(32)));

        gate.gateCandidates(selfOnly, selfRecord);

        // Self keeps a collapsed cigar and is never dropped.
        assertEquals("20M12S", selfOnly.get(0).LiftedCigar);
        assertFalse(selfOnly.get(0).Dropped);

        LiftedAlignment cleanSelf = selfAlignment(CHR_1, 100, "50M");
        LiftedAlignment collapsedAlt = refAlt(CHR_1, 300, "20M100N8M4S", true);
        List<LiftedAlignment> collapsedAltAlignments = new ArrayList<>(List.of(cleanSelf, collapsedAlt));
        SAMRecord collapseRecord = primaryRecord(CHR_1, 100, "50M");
        collapseRecord.setReadBases(bases("C".repeat(50)));

        gate.gateCandidates(collapsedAltAlignments, collapseRecord);

        // Collapsed XA alts are dropped from the alternate loci.
        assertTrue(collapsedAlt.Dropped);

    }

}
