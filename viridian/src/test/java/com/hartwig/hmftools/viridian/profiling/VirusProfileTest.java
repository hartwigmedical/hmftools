package com.hartwig.hmftools.viridian.profiling;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.viridian.profiling.common.GenomeStats;
import com.hartwig.hmftools.viridian.profiling.common.SummaryStats;
import com.hartwig.hmftools.viridian.profiling.select.RepresentativeGenome;
import com.hartwig.hmftools.viridian.profiling.select.RepresentativeStatus;
import com.hartwig.hmftools.viridian.reference.OncologyGroup;
import com.hartwig.hmftools.viridian.reference.VirusGenome;

import org.junit.Test;

public class VirusProfileTest
{
    private static final OncologyGroup GROUP_A = new OncologyGroup("Group A");
    private static final OncologyGroup GROUP_B = new OncologyGroup("Group B");
    private static final VirusGenome V1 = new VirusGenome("v1", 100, "Virus v1", GROUP_A);
    private static final VirusGenome W1 = new VirusGenome("w1", 100, "Virus w1", GROUP_B);

    @Test
    public void testFromPairsMeasurementWithStatus()
    {
        GenomeStats stats = stats(V1);

        List<VirusProfile> profiles = VirusProfile.from(List.of(representative(V1)), Map.of(V1, stats));

        assertEquals(List.of(new VirusProfile(stats, RepresentativeStatus.PASS)), profiles);
    }

    // Read assignment gave every read of V1 to another virus's representative, so V1 has no measurement and its oncology
    // group is dropped.
    @Test
    public void testFromRepresentativeWithoutReadsNotProfiled()
    {
        GenomeStats w1Stats = stats(W1);

        List<VirusProfile> profiles =
                VirusProfile.from(List.of(representative(V1), representative(W1)), Map.of(W1, w1Stats));

        assertEquals(List.of(new VirusProfile(w1Stats, RepresentativeStatus.PASS)), profiles);
    }

    private static RepresentativeGenome representative(VirusGenome genome)
    {
        return new RepresentativeGenome(genome, RepresentativeStatus.PASS);
    }

    private static GenomeStats stats(VirusGenome genome)
    {
        return new GenomeStats(
                genome, 10, 1, 80, SummaryStats.from(new int[genome.length()]), SummaryStats.from(new int[] { 100 }));
    }
}
