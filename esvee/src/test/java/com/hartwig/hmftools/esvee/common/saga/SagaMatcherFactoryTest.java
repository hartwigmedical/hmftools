package com.hartwig.hmftools.esvee.common.saga;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;

import static org.junit.Assert.assertEquals;

import java.util.Arrays;
import java.util.List;

import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.junit.Test;

import htsjdk.samtools.SAMSequenceDictionary;

public class SagaMatcherFactoryTest
{
    private static final int LOCATION_DISTANCE_MAX = 100;
    private static final String CHR1 = "chr1";
    private static final String CHR2 = "chr2";

    // subregion under test; expanded by LOCATION_DISTANCE_MAX → [900, 2100].
    private static final ChrBaseRegion SUBREGION = new ChrBaseRegion(CHR1, 1000, 2000);

    private static final SagaLocationMatcher.Config LOCATION_CONFIG =
            new SagaLocationMatcher.Config(LOCATION_DISTANCE_MAX);

    private int mVariantIdCounter = 0;

    // Creates a deletion assembly with breakend1 (FORWARD) at `position` on `chromosome`.
    // breakend2 is placed 10_000_000bp away so it never falls inside any test subregion.
    private SagaAssembly assemblyAt(final String chromosome, int position)
    {
        String id = "var" + (++mVariantIdCounter);
        SagaVariant variant = new SagaVariant(
                id,
                new SagaBreakend(new BasePosition(chromosome, position), FORWARD),
                new SagaBreakend(new BasePosition(chromosome, position + 10_000_000), REVERSE),
                ""
        );
        return new SagaAssembly(id, variant, List.of(50), "A".repeat(100));
    }

    private static SagaMatcherFactory factory(final SagaAssembly... assemblies)
    {
        SagaResource resource = new SagaResource("", Arrays.asList(assemblies), new SAMSequenceDictionary());
        return new SagaMatcherFactory(resource, LOCATION_CONFIG, SagaSequenceMatcherConfig.DEFAULT);
    }

    private static List<SagaVariant> slicedVariants(final SagaLocationMatcher matcher, final String chromosome)
    {
        List<SagaIndexedBreakend> breakends = matcher.getBreakends().getOrDefault(chromosome, List.of());
        return breakends.stream().map(SagaIndexedBreakend::variant).toList();
    }

    @Test
    public void testNoBreakendOnTargetChromosome()
    {
        // All breakends are on chr2; the subregion is on chr1 → no chr1 breakends in the matcher.
        SagaLocationMatcher matcher = factory(assemblyAt(CHR2, 1000), assemblyAt(CHR2, 1500))
                .createLocationMatcher(SUBREGION);

        assertEquals(List.of(), slicedVariants(matcher, CHR1));
    }

    @Test
    public void testNoBreakendWithinExpandedRegion()
    {
        // No breakend falls within the expanded region [900, 2100] → empty slice.
        SagaAssembly a100 = assemblyAt(CHR1, 100);
        SagaAssembly a200 = assemblyAt(CHR1, 200);
        SagaAssembly a5000 = assemblyAt(CHR1, 5000);
        SagaAssembly a6000 = assemblyAt(CHR1, 6000);

        SagaLocationMatcher matcher = factory(a100, a200, a5000, a6000).createLocationMatcher(SUBREGION);

        assertEquals(List.of(), slicedVariants(matcher, CHR1));
    }

    @Test
    public void testBreakendAtEdgeOfExpandedRegion()
    {
        // Breakends exactly on the expanded boundary (900 and 2100) are included.
        // Breakends just outside the boundary (899 and 2101) are excluded.
        SagaAssembly a899 = assemblyAt(CHR1, 899);
        SagaAssembly a900 = assemblyAt(CHR1, 900);
        SagaAssembly a2100 = assemblyAt(CHR1, 2100);
        SagaAssembly a2101 = assemblyAt(CHR1, 2101);

        SagaLocationMatcher matcher = factory(a899, a900, a2100, a2101).createLocationMatcher(SUBREGION);

        assertEquals(
                List.of(a900.variant(), a2100.variant()),
                slicedVariants(matcher, CHR1));
    }

    @Test
    public void testBreakendNearbyWithinTolerance()
    {
        // Breakends at 950 and 2050 lie outside the subregion [1000, 2000] but within
        // locationDistanceMax (100), so they fall inside the expanded region [900, 2100] → included.
        // Breakends outside the expanded region (100 and 5000) are excluded.
        SagaAssembly a100 = assemblyAt(CHR1, 100);
        SagaAssembly a950 = assemblyAt(CHR1, 950);
        SagaAssembly a2050 = assemblyAt(CHR1, 2050);
        SagaAssembly a5000 = assemblyAt(CHR1, 5000);

        SagaLocationMatcher matcher = factory(a100, a950, a2050, a5000).createLocationMatcher(SUBREGION);

        assertEquals(
                List.of(a950.variant(), a2050.variant()),
                slicedVariants(matcher, CHR1));
    }

    @Test
    public void testBreakendWithinSubregion()
    {
        // Breakends inside the subregion [1000, 2000] are included.
        // Breakends outside the expanded region [900, 2100] are excluded.
        SagaAssembly a100 = assemblyAt(CHR1, 100);
        SagaAssembly a1000 = assemblyAt(CHR1, 1000);
        SagaAssembly a1500 = assemblyAt(CHR1, 1500);
        SagaAssembly a2000 = assemblyAt(CHR1, 2000);
        SagaAssembly a5000 = assemblyAt(CHR1, 5000);

        SagaLocationMatcher matcher = factory(a100, a1000, a1500, a2000, a5000).createLocationMatcher(SUBREGION);

        assertEquals(
                List.of(a1000.variant(), a1500.variant(), a2000.variant()),
                slicedVariants(matcher, CHR1));
    }
}
