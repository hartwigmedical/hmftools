package com.hartwig.hmftools.compar;

import static com.hartwig.hmftools.common.driver.DriverCategory.ONCO;
import static com.hartwig.hmftools.common.driver.LikelihoodMethod.AMP;

import static junit.framework.TestCase.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.common.driver.DriverType;
import com.hartwig.hmftools.common.driver.ImmutableDriverCatalog;
import com.hartwig.hmftools.compar.common.CommonUtils;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.MismatchType;
import com.hartwig.hmftools.compar.driver.DriverComparer;
import com.hartwig.hmftools.compar.driver.DriverData;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class DriverTest
{
    @Test
    public void testDriverDiffsWithoutMatches()
    {
        List<Mismatch> mismatches = generateTestMismatches(false);

        assertEquals(3, mismatches.size());
        assertEquals(1, mismatches.stream().filter(x -> x.MismatchType == MismatchType.REF_ONLY).count());
        assertEquals(1, mismatches.stream().filter(x -> x.MismatchType == MismatchType.NEW_ONLY).count());
        assertEquals(1, mismatches.stream().filter(x -> x.MismatchType == MismatchType.VALUE).count());
    }

    @Test
    public void testDriverDiffsWithMatches()
    {
        List<Mismatch> mismatches = generateTestMismatches(true);

        assertEquals(4, mismatches.size());
        assertEquals(1, mismatches.stream().filter(x -> x.MismatchType == MismatchType.REF_ONLY).count());
        assertEquals(1, mismatches.stream().filter(x -> x.MismatchType == MismatchType.NEW_ONLY).count());
        assertEquals(1, mismatches.stream().filter(x -> x.MismatchType == MismatchType.VALUE).count());
        assertEquals(1, mismatches.stream().filter(x -> x.MismatchType == MismatchType.FULL_MATCH).count());
    }

    @NotNull
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
