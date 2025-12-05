package com.hartwig.hmftools.compar.driver;

import static com.hartwig.hmftools.common.driver.DriverCategory.ONCO;
import static com.hartwig.hmftools.common.driver.LikelihoodMethod.AMP;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.compar.common.CommonUtils.FLD_CHROMOSOME_BAND;
import static com.hartwig.hmftools.compar.driver.DriverData.FLD_LIKELIHOOD;
import static com.hartwig.hmftools.compar.driver.DriverData.FLD_LIKE_METHOD;
import static com.hartwig.hmftools.compar.driver.DriverData.FLD_MAX_COPY_NUMBER;
import static com.hartwig.hmftools.compar.driver.DriverData.FLD_MIN_COPY_NUMBER;

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
        comparer = new DriverComparer(new ComparConfig());
        builder = TestDriverDataBuilder.BUILDER;
        DriverData alternateValueSource = builder.createWithAlternateDefaults();
        fieldToAlternateValueInitializer = Map.of(
                FLD_LIKE_METHOD, b -> b.likelihoodMethod = alternateValueSource.DriverCatalog.likelihoodMethod(),
                FLD_LIKELIHOOD, b -> b.likelihood = alternateValueSource.DriverCatalog.driverLikelihood(),
                FLD_MIN_COPY_NUMBER, b -> b.minCopyNumber = alternateValueSource.DriverCatalog.minCopyNumber(),
                FLD_MAX_COPY_NUMBER, b -> b.maxCopyNumber = alternateValueSource.DriverCatalog.maxCopyNumber(),
                FLD_CHROMOSOME, b -> {
                    b.chromosome = alternateValueSource.DriverCatalog.chromosome();
                    b.comparisonChromosome = alternateValueSource.mComparisonChromosome;
                },
                FLD_CHROMOSOME_BAND, b -> b.chromosomeBand = alternateValueSource.DriverCatalog.chromosomeBand()
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
        assertEquals(1, mismatches.stream().filter(x -> x.MismatchType() == MismatchType.REF_ONLY).count());
        assertEquals(1, mismatches.stream().filter(x -> x.MismatchType() == MismatchType.NEW_ONLY).count());
        assertEquals(1, mismatches.stream().filter(x -> x.MismatchType() == MismatchType.VALUE).count());
    }

    @Test
    public void testDriverDiffsWithMatches()
    {
        List<Mismatch> mismatches = generateTestMismatches(true);

        assertEquals(4, mismatches.size());
        assertEquals(1, mismatches.stream().filter(x -> x.MismatchType() == MismatchType.REF_ONLY).count());
        assertEquals(1, mismatches.stream().filter(x -> x.MismatchType() == MismatchType.NEW_ONLY).count());
        assertEquals(1, mismatches.stream().filter(x -> x.MismatchType() == MismatchType.VALUE).count());
        assertEquals(1, mismatches.stream().filter(x -> x.MismatchType() == MismatchType.FULL_MATCH).count());
    }

    private static List<Mismatch> generateTestMismatches(final boolean includeMatches)
    {
        ComparConfig config = new ComparConfig();
        DriverComparer driverComparer = new DriverComparer(config);

        driverComparer.registerThresholds(config.Thresholds);

        List<ComparableItem> refItems = Lists.newArrayList();
        List<ComparableItem> newItems = Lists.newArrayList();

        refItems.add(new DriverData(createDriverCatalog("AR", DriverType.AMP, 1.0, 6), "1", false));

        newItems.add(new DriverData(createDriverCatalog("TP53", DriverType.DEL, 1.0, 0.2), "2", false));

        refItems.add(new DriverData(createDriverCatalog("KRAS", DriverType.MUTATION, 0.7, 2), "3", false));
        newItems.add(new DriverData(createDriverCatalog("KRAS", DriverType.MUTATION, 0.5, 2), "3", false));

        refItems.add(new DriverData(createDriverCatalog("BRAF", DriverType.HOM_DEL_DISRUPTION, 0.9, 2), "4", false));
        newItems.add(new DriverData(createDriverCatalog("BRAF", DriverType.HOM_DEL_DISRUPTION, 0.9, 2), "4", false));

        List<Mismatch> mismatches = Lists.newArrayList();
        CommonUtils.compareItems(mismatches, MatchLevel.REPORTABLE, config.Thresholds, includeMatches, refItems, newItems);
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
