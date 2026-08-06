package com.hartwig.hmftools.compar.driver;

import static com.hartwig.hmftools.common.driver.DriverCategory.ONCO;
import static com.hartwig.hmftools.common.driver.LikelihoodMethod.AMP;
import static com.hartwig.hmftools.compar.driver.TestDriverDataBuilder.buildPurityData;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import static junit.framework.TestCase.assertEquals;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.common.driver.DriverType;
import com.hartwig.hmftools.common.driver.ImmutableDriverCatalog;
import com.hartwig.hmftools.common.purple.PurplePurity;
import com.hartwig.hmftools.common.purple.ReportedStatus;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.ComparableItemTest;
import com.hartwig.hmftools.compar.common.CommonUtils;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.MismatchType;

import org.junit.Before;
import org.junit.Test;

public class DriverDataTest extends ComparableItemTest<DriverData, DriverComparer, TestDriverDataBuilder>
{
    @Before
    public void setUp()
    {
        comparer = new DriverComparer(new ComparConfig(), Collections.emptyMap());
        builder = TestDriverDataBuilder.BUILDER;
        DriverData alternateValueSource = builder.createWithAlternateDefaults();
        fieldToAlternateValueInitializer = Map.of(
                DriverComparer.Fields.LikelihoodMethod.toString(), b -> b.likelihoodMethod = alternateValueSource.DriverCatalog.likelihoodMethod(),
                DriverComparer.Fields.Likelihood.toString(), b -> b.likelihood = alternateValueSource.DriverCatalog.driverLikelihood(),
                DriverComparer.Fields.MinCopyNumber.toString(), b -> b.minCopyNumber = alternateValueSource.DriverCatalog.minCopyNumber(),
                DriverComparer.Fields.MaxCopyNumber.toString(), b -> b.maxCopyNumber = alternateValueSource.DriverCatalog.maxCopyNumber(),
                DriverComparer.Fields.Chromosome.toString(), b -> b.comparisonChromosome = alternateValueSource.ComparisonChromosome,
                DriverComparer.Fields.ChromosomeBand.toString(), b -> b.chromosomeBand = alternateValueSource.DriverCatalog.chromosomeBand()
        );
        nameToAlternateIndexInitializer = Map.of(
                "gene", b -> b.gene = alternateValueSource.DriverCatalog.gene(),
                "driver", b -> b.driver = alternateValueSource.DriverCatalog.driver(),
                "nonCanonicalTranscript", b -> {
                    b.transcript = alternateValueSource.DriverCatalog.transcript();
                    b.isCanonical = alternateValueSource.DriverCatalog.isCanonical();
                }
        );
        reportabilityFieldToFalseReportabilityInitializer = Collections.emptyMap();
        nameToNonPassInitializer = Map.of("nonPass", b -> b.isPass = false);
    }
    
    @Test
    public void checkTranscriptWhenNeeded()
    {
        DriverData victim = TestDriverDataBuilder.BUILDER.create();
        DriverData victimOtherTranscript = TestDriverDataBuilder.BUILDER.create(b -> b.transcript = "OTHER");
        DriverData victimCanonical = TestDriverDataBuilder.BUILDER.create(b -> b.isCanonical = true);
        DriverData victimOtherTranscriptAndCanonical = TestDriverDataBuilder.BUILDER.create(b -> {
            b.transcript = "OTHER";
            b.isCanonical = true;
        });
        DriverData victimWithoutCheckTranscript = TestDriverDataBuilder.BUILDER.create(b ->
        {
            b.transcript = "OTHER";
            b.isCanonical = true;
            b.checkTranscript = false;
        });

        assertTrue(victim.matches(victim));

        assertTrue(victim.matches(victimOtherTranscript));
        assertTrue(victim.matches(victimCanonical));
        assertTrue(victimOtherTranscript.matches(victim));
        assertTrue(victimCanonical.matches(victim));

        assertFalse(victim.matches(victimOtherTranscriptAndCanonical));

        assertTrue(victimOtherTranscriptAndCanonical.matches(victimWithoutCheckTranscript));
        assertTrue(victimWithoutCheckTranscript.matches(victimOtherTranscriptAndCanonical));

        assertFalse(victim.matches(victimWithoutCheckTranscript));
        assertFalse(victimWithoutCheckTranscript.matches(victim));
    }

    @Test
    public void doNotCheckTranscriptWhenNotNeeded()
    {
        DriverData victim = TestDriverDataBuilder.BUILDER.create(b -> b.checkTranscript = false);
        DriverData victimOtherTranscript = TestDriverDataBuilder.BUILDER.create(b ->
        {
            b.transcript = "OTHER";
            b.checkTranscript = false;
        });
        DriverData victimCanonical = TestDriverDataBuilder.BUILDER.create(b ->
        {
            b.isCanonical = true;
            b.checkTranscript = false;
        });
        DriverData victimOtherTranscriptAndCanonical = TestDriverDataBuilder.BUILDER.create(b ->
        {
            b.transcript = "OTHER";
            b.isCanonical = true;
            b.checkTranscript = false;
        });

        assertTrue(victim.matches(victim));
        assertTrue(victim.matches(victimOtherTranscript));
        assertTrue(victim.matches(victimCanonical));
        assertTrue(victim.matches(victimOtherTranscriptAndCanonical));
        assertTrue(victimOtherTranscript.matches(victim));
        assertTrue(victimCanonical.matches(victim));
        assertTrue(victimOtherTranscriptAndCanonical.matches(victim));
    }

    @Test
    public void testDriverDiffsWithoutMatches()
    {
        List<Mismatch> mismatches = generateTestMismatches(false);

        assertEquals(3, mismatches.size());
        assertEquals(1, mismatches.stream().filter(x -> x.Type == MismatchType.OLD_ONLY).count());
        assertEquals(1, mismatches.stream().filter(x -> x.Type == MismatchType.NEW_ONLY).count());
        assertEquals(1, mismatches.stream().filter(x -> x.Type == MismatchType.VALUE).count());
    }

    @Test
    public void testDriverDiffsWithMatches()
    {
        List<Mismatch> mismatches = generateTestMismatches(true);

        assertEquals(4, mismatches.size());
        assertEquals(1, mismatches.stream().filter(x -> x.Type == MismatchType.OLD_ONLY).count());
        assertEquals(1, mismatches.stream().filter(x -> x.Type == MismatchType.NEW_ONLY).count());
        assertEquals(1, mismatches.stream().filter(x -> x.Type == MismatchType.VALUE).count());
        assertEquals(1, mismatches.stream().filter(x -> x.Type == MismatchType.FULL_MATCH).count());
    }

    private static List<Mismatch> generateTestMismatches(final boolean includeMatches)
    {
        ComparConfig config = new ComparConfig();
        DriverComparer driverComparer = new DriverComparer(config, Collections.emptyMap());

        List<ComparableItem> refItems = Lists.newArrayList();
        List<ComparableItem> newItems = Lists.newArrayList();

        PurplePurity purity = buildPurityData();

        refItems.add(
                new DriverData(createDriverCatalog("AR", DriverType.AMP, 1.0, 6),
                        purity, "1", false, true, driverComparer.fieldsList()));

        newItems.add(
                new DriverData(createDriverCatalog("TP53", DriverType.DEL, 1.0, 0.2),
                        purity, "2", false, true, driverComparer.fieldsList()));

        refItems.add(new DriverData(createDriverCatalog("KRAS", DriverType.MUTATION, 0.7, 2),
                purity, "3", false, true, driverComparer.fieldsList()));
        newItems.add(new DriverData(createDriverCatalog("KRAS", DriverType.MUTATION, 0.5, 2),
                purity, "3", false, true, driverComparer.fieldsList()));

        refItems.add(new DriverData(createDriverCatalog("BRAF", DriverType.HOM_DEL_DISRUPTION, 0.9, 2),
                purity, "4", false, true, driverComparer.fieldsList()));

        newItems.add(new DriverData(createDriverCatalog("BRAF", DriverType.HOM_DEL_DISRUPTION, 0.9, 2),
                purity, "4", false, true, driverComparer.fieldsList()));

        List<Mismatch> mismatches = Lists.newArrayList();
        CommonUtils.compareItems(driverComparer, mismatches, MatchLevel.REPORTABLE, includeMatches, refItems, newItems);
        return mismatches;
    }

    private static DriverCatalog createDriverCatalog(final String gene, final DriverType type, double likelihood, double minCopyNumber)
    {
        return ImmutableDriverCatalog.builder()
                .chromosome("1")
                .chromosomeBand("q28")
                .gene(gene)
                .transcript("")
                .isCanonical(true)
                .driver(type)
                .category(ONCO)
                .likelihoodMethod(AMP)
                .reportedStatus(ReportedStatus.REPORTED)
                .driverLikelihood(likelihood)
                .missense(0)
                .nonsense(0)
                .splice(0)
                .inframe(0)
                .frameshift(0)
                .biallelic(false)
                .minCopyNumber(minCopyNumber)
                .maxCopyNumber(minCopyNumber).build();
    }
}
