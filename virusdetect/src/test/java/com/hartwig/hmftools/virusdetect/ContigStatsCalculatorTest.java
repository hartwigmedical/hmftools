package com.hartwig.hmftools.virusdetect;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Map;

import org.junit.Test;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

public class ContigStatsCalculatorTest
{
    private static final double EPSILON = 1e-9;

    // v1 (length 20): read r1 covers 1-10 (score 10), read r2 covers 6-10 (score 8); a lower-scoring second alignment of
    // r1 on v1 is deduped for depth/score but still counts toward r1's two alignments there. v2 (length 10): r3 full.
    @Test
    public void testComputesPerContigDepthAndCoverage()
    {
        ViralReference reference = reference();

        List<ReadAlignment> alignments = List.of(
                alignment("r1", "v1", 1, 10, 10, 0),
                alignment("r1", "v1", 1, 10, 5, 0),   // second, lower-scoring alignment of r1 on v1
                alignment("r2", "v1", 6, 5, 8, 0),
                alignment("r3", "v2", 1, 10, 9, 0));

        Map<String, ContigStats> stats = new ContigStatsCalculator().compute(alignments, reference);

        assertEquals(2, stats.size());

        ContigStats v1 = stats.get("v1");
        assertEquals(20, v1.contigLength());
        assertEquals(2, v1.readCount());              // r1 counted once despite two alignments
        assertEquals(1, v1.multiAlignReads());          // r1 has two alignments on v1, r2 one
        assertEquals(1.5, v1.alignPerRead().mean(), EPSILON);   // r1 -> 2, r2 -> 1
        assertEquals(1.0, v1.alignPerRead().min(), EPSILON);
        assertEquals(2.0, v1.alignPerRead().max(), EPSILON);
        assertEquals(10, v1.coveredBases());          // positions 1-10
        assertEquals(0.0, v1.depth().min(), EPSILON);        // positions 11-20 uncovered
        assertEquals(2.0, v1.depth().max(), EPSILON);        // positions 6-10 covered by both reads
        assertEquals(15.0 / 20, v1.depth().mean(), EPSILON);
        assertEquals(0.0, v1.depth().p50(), EPSILON);        // >half the contig uncovered, so median depth is 0
        assertEquals(2.0, v1.depth().p95(), EPSILON);
        assertEquals(9.0, v1.alignerScore().mean(), EPSILON);   // (10 + 8) / 2
        assertEquals(8.0, v1.alignerScore().min(), EPSILON);
        assertEquals(10.0, v1.alignerScore().max(), EPSILON);
        assertEquals(0.5, v1.coverageFraction(), EPSILON);

        ContigStats v2 = stats.get("v2");
        assertEquals(1, v2.readCount());
        assertEquals(0, v2.multiAlignReads());          // r3 has a single alignment on v2
        assertEquals(1.0, v2.alignPerRead().mean(), EPSILON);
        assertEquals(10, v2.coveredBases());
        assertEquals(1.0, v2.depth().min(), EPSILON);
        assertEquals(1.0, v2.depth().max(), EPSILON);
        assertEquals(1.0, v2.depth().mean(), EPSILON);
        assertEquals(9.0, v2.alignerScore().mean(), EPSILON);
        assertEquals(1.0, v2.coverageFraction(), EPSILON);
    }

    // r1 and r2 each align to both contigs, more closely to v1 (lower divergence). Both reads' votes and their contested
    // margins land on v1; v2, never a read's best, holds no margins.
    @Test
    public void testVotesAndMarginsAttributeStrainSupport()
    {
        ViralReference reference = reference();

        List<ReadAlignment> alignments = List.of(
                alignment("r1", "v1", 1, 10, 10, 1),
                alignment("r1", "v2", 1, 10, 6, 4),
                alignment("r2", "v1", 1, 10, 10, 0),
                alignment("r2", "v2", 1, 10, 5, 5));

        // Injected correct-base probability 0.5, so each extra divergent base halves a contig's weight (0.5^diff).
        Map<String, ContigStats> stats = new ContigStatsCalculator(0.5).compute(alignments, reference);

        ContigStats v1 = stats.get("v1");
        ContigStats v2 = stats.get("v2");

        // r1: v1 weight 0.5^0, v2 0.5^3; r2: v1 0.5^0, v2 0.5^5. Votes per read sum to 1, so both contigs sum to 2.
        assertEquals(1.0 / (1 + Math.pow(0.5, 3)) + 1.0 / (1 + Math.pow(0.5, 5)), v1.readVotes(), EPSILON);
        assertEquals(2.0 - v1.readVotes(), v2.readVotes(), EPSILON);

        // v1 wins both reads; margins are runner-up minus best divergence: r1 4-1=3, r2 5-0=5.
        assertEquals(2, v1.readsBestInRivals());
        SummaryStats v1Margins = v1.margins().orElseThrow();
        assertEquals(4.0, v1Margins.mean(), EPSILON);
        assertEquals(3.0, v1Margins.min(), EPSILON);
        assertEquals(3.0, v1Margins.p50(), EPSILON);
        assertEquals(5.0, v1Margins.p95(), EPSILON);
        assertEquals(5.0, v1Margins.max(), EPSILON);

        // v2 is never a read's best, so it holds no margins.
        assertEquals(0, v2.readsBestInRivals());
        assertTrue(v2.margins().isEmpty());
    }

    // A read that ties across contigs has no strict winner, so no contig is credited with a best-in-rivals read; the
    // tied contigs stay symmetric. The read still splits its vote evenly between them.
    @Test
    public void testTiedReadCreditsNoContigAsBest()
    {
        ViralReference reference = reference();

        List<ReadAlignment> alignments = List.of(
                alignment("r", "v1", 1, 10, 10, 2),
                alignment("r", "v2", 1, 10, 10, 2));

        Map<String, ContigStats> stats = new ContigStatsCalculator(0.5).compute(alignments, reference);

        assertEquals(0, stats.get("v1").readsBestInRivals());
        assertEquals(0, stats.get("v2").readsBestInRivals());
        assertTrue(stats.get("v1").margins().isEmpty());
        assertTrue(stats.get("v2").margins().isEmpty());
        assertEquals(0.5, stats.get("v1").readVotes(), EPSILON);
        assertEquals(0.5, stats.get("v2").readVotes(), EPSILON);
    }

    // Alignments whose clip projects past a contig end straddle the circular genome origin: they are dropped from the
    // stats and counted separately. A clip that stays within the contig is kept.
    @Test
    public void testDropsAlignmentsClippingOverContigEnds()
    {
        ViralReference reference = reference();   // v1 length 20

        List<ReadAlignment> alignments = List.of(
                clipped("r1", "v1", 1, 10, 10, 0),    // left clip projects to -9: over the start
                clipped("r2", "v1", 11, 20, 0, 10),   // right clip projects to 30: over the end
                alignment("r3", "v1", 6, 10, 9, 0),   // no clip: kept
                clipped("r4", "v1", 10, 19, 5, 0));   // clip projects to 5, within the contig: kept

        Map<String, ContigStats> stats = new ContigStatsCalculator().compute(alignments, reference);

        ContigStats v1 = stats.get("v1");
        assertEquals(2, v1.readCount());            // r3 and r4 kept
        assertEquals(2, v1.originClippedReads());   // r1 and r2 dropped
    }

    // A clean, full-length alignment covering [start, start + length - 1] with no clips.
    private static ReadAlignment alignment(String readName, String contig, int start, int length, int alignerScore, int divergence)
    {
        return new ReadAlignment(readName, contig, start, start + length - 1, 0, 0, alignerScore, divergence,
                List.of(new ReadAlignment.AlignedInterval(start, length)));
    }

    // An alignment covering [start, end] with the given clip lengths hanging off each side.
    private static ReadAlignment clipped(String readName, String contig, int start, int end, int leftClip, int rightClip)
    {
        return new ReadAlignment(readName, contig, start, end, leftClip, rightClip, 10, leftClip + rightClip,
                List.of(new ReadAlignment.AlignedInterval(start, end - start + 1)));
    }

    private static ViralReference reference()
    {
        List<ViralContig> contigs = List.of(
                new ViralContig("v1", 20, "Virus 1", "Group 1"),
                new ViralContig("v2", 10, "Virus 2", "Group 2"));
        SAMSequenceDictionary dictionary = new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("v1", 20),
                new SAMSequenceRecord("v2", 10)));
        return new ViralReference(contigs, dictionary);
    }
}
