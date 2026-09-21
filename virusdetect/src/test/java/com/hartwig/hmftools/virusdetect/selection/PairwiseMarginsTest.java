package com.hartwig.hmftools.virusdetect.selection;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.virusdetect.detection.read_align.AlignedInterval;
import com.hartwig.hmftools.virusdetect.detection.read_align.ViralReadAlignment;
import com.hartwig.hmftools.virusdetect.detection.read_align.ViralReadAlignments;
import com.hartwig.hmftools.virusdetect.reference.OncologyGroup;
import com.hartwig.hmftools.virusdetect.reference.ViralContig;

import org.junit.Test;

public class PairwiseMarginsTest
{
    private static final int LENGTH = 100;
    private static final double MEAN_READ_LENGTH = 150.0;

    private static final OncologyGroup GROUP_A = new OncologyGroup("Group A");
    private static final OncologyGroup GROUP_H = new OncologyGroup("Group H");

    private static final ViralContig V1 = new ViralContig("v1", LENGTH, "Virus 1", GROUP_A);
    private static final ViralContig V2 = new ViralContig("v2", LENGTH, "Virus 2", GROUP_A);
    private static final ViralContig H1 = new ViralContig("h1", LENGTH, "Virus H", GROUP_H);

    @Test
    public void testPairsContigsWithinGroupByWinningMargin()
    {
        List<ViralReadAlignment> alignments = List.of(
                alignment("r1", V1, 2),
                alignment("r1", V2, 7),   // v1 wins r1 by 5
                alignment("r2", V1, 0),
                alignment("r2", V2, 3),   // v1 wins r2 by 3
                alignment("r3", V1, 1),
                alignment("r3", H1, 1));  // cross-group: not paired

        PairwiseMargins margins = PairwiseMargins.from(ViralReadAlignments.from(alignments, MEAN_READ_LENGTH));

        // Both reads shared between v1 and v2, both directions
        assertEquals(2, margins.sharedReads(V1, V2));
        assertEquals(2, margins.sharedReads(V2, V1));

        // v1's winning margins over v2 are 5 (r1) and 3 (r2)
        assertEquals(2, margins.readsWinningBy(V1, V2, 3));
        assertEquals(1, margins.readsWinningBy(V1, V2, 5));
        assertEquals(0, margins.readsWinningBy(V1, V2, 6));

        // v2 never fits a shared read better than v1
        assertEquals(0, margins.readsWinningBy(V2, V1, 1));

        // The cross-group read produced no pair
        assertEquals(0, margins.sharedReads(V1, H1));
    }

    private static ViralReadAlignment alignment(String readName, ViralContig contig, int divergence)
    {
        return new ViralReadAlignment(
                readName, contig, 1, LENGTH, 0, 0, 100, divergence, List.of(new AlignedInterval(1, LENGTH)));
    }
}
