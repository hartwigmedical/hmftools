package com.hartwig.hmftools.virusdetect;

import static org.junit.Assert.assertEquals;

import java.util.List;

import org.junit.Test;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

public class PairwiseMarginsTest
{
    private static final ViralReference REFERENCE = reference();

    // v1 and v2 share a group; h1 is in another. Two reads align to both v1 and v2, each fitting v1 better; a third
    // spans the group boundary and must not be paired.
    @Test
    public void testPairsContigsWithinGroupByWinningMargin()
    {
        ViralContig v1 = REFERENCE.contig("v1");
        ViralContig v2 = REFERENCE.contig("v2");

        List<ViralAlignment> alignments = List.of(
                alignment("r1", "v1", 2),
                alignment("r1", "v2", 7),   // v1 wins r1 by 5
                alignment("r2", "v1", 0),
                alignment("r2", "v2", 3),   // v1 wins r2 by 3
                alignment("r3", "v1", 1),
                alignment("r3", "h1", 1));  // cross-group: not paired

        PairwiseMargins margins = PairwiseMargins.from(ViralAlignments.from(alignments, 150.0));

        // Both reads shared between v1 and v2, both directions
        assertEquals(2, margins.sharedReads(v1, v2));
        assertEquals(2, margins.sharedReads(v2, v1));

        // v1's winning margins over v2 are 5 (r1) and 3 (r2)
        assertEquals(2, margins.readsWinningBy(v1, v2, 3));
        assertEquals(1, margins.readsWinningBy(v1, v2, 5));
        assertEquals(0, margins.readsWinningBy(v1, v2, 6));

        // v2 never fits a shared read better than v1
        assertEquals(0, margins.readsWinningBy(v2, v1, 1));

        // The cross-group read produced no pair
        assertEquals(0, margins.sharedReads(v1, REFERENCE.contig("h1")));
    }

    // An alignment whose clip projects past a contig end straddles the circular origin and is dropped before pairing.
    @Test
    public void testDropsOriginStraddlers()
    {
        List<ViralAlignment> alignments = List.of(
                alignment("r1", "v1", 2),
                straddler("r1", "v2"));   // dropped, so r1 has only v1 left and forms no pair

        PairwiseMargins margins = PairwiseMargins.from(ViralAlignments.from(alignments, 150.0));

        assertEquals(0, margins.sharedReads(REFERENCE.contig("v1"), REFERENCE.contig("v2")));
    }

    private static ViralAlignment alignment(String readName, String contig, int divergence)
    {
        return new ViralAlignment(
                readName, REFERENCE.contig(contig), 1, 100, 0, 0, 100, divergence,
                List.of(new ViralAlignment.AlignedInterval(1, 100)));
    }

    // A right clip projecting well past the 100-base contig end.
    private static ViralAlignment straddler(String readName, String contig)
    {
        return new ViralAlignment(
                readName, REFERENCE.contig(contig), 1, 100, 0, 50, 100, 50,
                List.of(new ViralAlignment.AlignedInterval(1, 100)));
    }

    private static ViralReference reference()
    {
        List<ViralContig> contigs = List.of(
                new ViralContig("v1", 100, "Virus 1", "Group A"),
                new ViralContig("v2", 100, "Virus 2", "Group A"),
                new ViralContig("h1", 100, "Virus H", "Group H"));
        SAMSequenceDictionary dictionary = new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("v1", 100),
                new SAMSequenceRecord("v2", 100),
                new SAMSequenceRecord("h1", 100)));
        return new ViralReference(contigs, dictionary);
    }
}
