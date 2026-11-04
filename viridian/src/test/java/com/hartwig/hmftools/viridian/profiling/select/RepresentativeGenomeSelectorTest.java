package com.hartwig.hmftools.viridian.profiling.select;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.viridian.profiling.common.GenomeStats;
import com.hartwig.hmftools.viridian.profiling.common.SummaryStats;
import com.hartwig.hmftools.viridian.profiling.support.GenomeSupport;
import com.hartwig.hmftools.viridian.reference.Virus;
import com.hartwig.hmftools.viridian.reference.VirusGenome;
import com.hartwig.hmftools.viridian.reference.VirusGroup;

import org.junit.Test;

public class RepresentativeGenomeSelectorTest
{
    private static final int LENGTH = 1000;

    private static final VirusGroup GROUP_A = new VirusGroup("Group A");
    private static final VirusGroup GROUP_B = new VirusGroup("Group B");
    private static final VirusGroup GROUP_H = new VirusGroup("Group H");

    private static final VirusGenome V1 = new VirusGenome("v1", LENGTH, new Virus("Virus v1", GROUP_A));
    private static final VirusGenome V2 = new VirusGenome("v2", LENGTH, new Virus("Virus v2", GROUP_A));
    private static final VirusGenome V3 = new VirusGenome("v3", LENGTH, new Virus("Virus v3", GROUP_A));
    private static final VirusGenome V4 = new VirusGenome("v4", LENGTH, new Virus("Virus v4", GROUP_A));
    private static final VirusGenome W1 = new VirusGenome("w1", LENGTH, new Virus("Virus w1", GROUP_B));
    private static final VirusGenome H1 = new VirusGenome("h1", LENGTH, new Virus("Virus h1", GROUP_H));

    @Test
    public void testSelectRanksGenomesPerVirusGroup()
    {
        // V3 ranks 4th in group A, below the secondary genome count. W1 has few votes, but any read makes its group present.
        List<GenomeSupport> genomeSupports = List.of(
                withReads(V1, 100), withReads(V2, 120), withReads(V3, 5), withReads(V4, 50), withReads(W1, 1));

        List<RepresentativeSelection> selections = RepresentativeGenomeSelector.select(genomeSupports);

        assertEquals(
                List.of(new RepresentativeSelection(V2, List.of(V1, V4)), new RepresentativeSelection(W1, List.of())),
                selections);
    }

    @Test
    public void testSelectSkipsGenomeWithNoReads()
    {
        // H1's only reads clipped over its origin and were dropped.
        SummaryStats noDepth = SummaryStats.from(new int[LENGTH]);
        GenomeSupport noReads = new GenomeSupport(new GenomeStats(H1, 0, 2, 0, noDepth, null), 0, null, 0.0, 0);

        List<RepresentativeSelection> selections = RepresentativeGenomeSelector.select(List.of(withReads(W1, 10), noReads));

        assertEquals(List.of(new RepresentativeSelection(W1, List.of())), selections);
    }

    @Test
    public void testSelectTiedVotesByGenomeName()
    {
        List<RepresentativeSelection> selections =
                RepresentativeGenomeSelector.select(List.of(withReads(V2, 100), withReads(V1, 100)));

        assertEquals(List.of(new RepresentativeSelection(V1, List.of(V2))), selections);
    }

    // Selection reads only the virus genome, its read count and its read votes. The remaining fields are filler,
    // set to values consistent with a present virus genome so nothing reads as contradictory.
    private static GenomeSupport withReads(VirusGenome genome, double votes)
    {
        SummaryStats filler = SummaryStats.from(new int[] { 1 });
        return new GenomeSupport(new GenomeStats(genome, 100, 0, LENGTH / 2, filler, filler), 0, filler, votes, 1);
    }
}
