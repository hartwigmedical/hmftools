package com.hartwig.hmftools.viridian.profiling.select;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.viridian.profiling.common.GenomeStats;
import com.hartwig.hmftools.viridian.profiling.common.SummaryStats;
import com.hartwig.hmftools.viridian.profiling.support.GenomeFilterStatus;
import com.hartwig.hmftools.viridian.profiling.support.GenomeSupport;
import com.hartwig.hmftools.viridian.reference.OncologyGroup;
import com.hartwig.hmftools.viridian.reference.VirusGenome;

import org.junit.Test;

public class RepresentativeGenomeSelectorTest
{
    private static final int LENGTH = 1000;

    private static final OncologyGroup GROUP_A = new OncologyGroup("Group A");
    private static final OncologyGroup GROUP_B = new OncologyGroup("Group B");
    private static final OncologyGroup GROUP_H = new OncologyGroup("Group H");

    private static final VirusGenome V1 = new VirusGenome("v1", LENGTH, "Virus v1", GROUP_A);
    private static final VirusGenome V2 = new VirusGenome("v2", LENGTH, "Virus v2", GROUP_A);
    private static final VirusGenome V3 = new VirusGenome("v3", LENGTH, "Virus v3", GROUP_A);
    private static final VirusGenome W1 = new VirusGenome("w1", LENGTH, "Virus w1", GROUP_B);
    private static final VirusGenome H1 = new VirusGenome("h1", LENGTH, "Virus h1", GROUP_H);

    @Test
    public void testSelectTopVotesGenomePerOncologyGroup()
    {
        // V3 has the most votes in group A but was filtered; group H has only a filtered genome.
        List<GenomeSupport> genomeSupports = List.of(
                passed(V1, 100), passed(V2, 120), filtered(V3, 500), passed(W1, 10), filtered(H1, 100));

        List<VirusGenome> representatives = RepresentativeGenomeSelector.select(genomeSupports);

        assertEquals(List.of(V2, W1), representatives);
    }

    @Test
    public void testSelectTiedVotesByGenomeName()
    {
        List<VirusGenome> representatives =
                RepresentativeGenomeSelector.select(List.of(passed(V2, 100), passed(V1, 100)));

        assertEquals(List.of(V1), representatives);
    }

    private static GenomeSupport passed(VirusGenome genome, double votes)
    {
        return genomeSupport(genome, GenomeFilterStatus.PASS, votes);
    }

    private static GenomeSupport filtered(VirusGenome genome, double votes)
    {
        return genomeSupport(genome, GenomeFilterStatus.LOW_COVERAGE, votes);
    }

    // Selection reads only the virus genome, its filter status and its read votes. The remaining fields are filler,
    // set to values consistent with a present virus genome so nothing reads as contradictory.
    private static GenomeSupport genomeSupport(VirusGenome genome, GenomeFilterStatus filterStatus, double votes)
    {
        SummaryStats filler = SummaryStats.from(new int[] { 1 });
        return new GenomeSupport(
                new GenomeStats(genome, 100, 0, LENGTH / 2, filler, filler), filterStatus, 0, filler, votes);
    }
}
