package com.hartwig.hmftools.compar;

import static com.hartwig.hmftools.common.drivercatalog.DriverCategory.ONCO;
import static com.hartwig.hmftools.common.drivercatalog.DriverType.checkConvertType;
import static com.hartwig.hmftools.common.drivercatalog.LikelihoodMethod.AMP;

import static junit.framework.TestCase.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.common.drivercatalog.ImmutableDriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.LikelihoodMethod;
import com.hartwig.hmftools.compar.driver.DriverComparer;
import com.hartwig.hmftools.compar.driver.DriverData;

import org.junit.Test;

public class DriverTest
{
    @Test
    public void testDriverDiffs()
    {

        List<Mismatch> mismatches = Lists.newArrayList();

        ComparConfig config = new ComparConfig();
        DriverComparer driverComparer = new DriverComparer(config);

        driverComparer.registerThresholds(config.Thresholds);

        List<ComparableItem> refItems = Lists.newArrayList();
        List<ComparableItem> newItems = Lists.newArrayList();

        refItems.add(new DriverData(createDriverCatalog("AR", DriverType.AMP, 1.0, 6)));

        newItems.add(new DriverData(createDriverCatalog("TP53", DriverType.DEL, 1.0, 0.2)));

        refItems.add(new DriverData(createDriverCatalog("KRAS", DriverType.MUTATION, 0.7, 2)));
        newItems.add(new DriverData(createDriverCatalog("KRAS", DriverType.MUTATION, 0.5, 2)));

        refItems.add(new DriverData(createDriverCatalog("CDKN2A", DriverType.MUTATION, 1.0, 2)));
        newItems.add(new DriverData(createDriverCatalog("CDKN2A", DriverType.MUTATION, 1.0, 1.5)));

        CommonUtils.compareItems(mismatches, MatchLevel.REPORTABLE, config.Thresholds, refItems, newItems);

        assertEquals(4, mismatches.size());
        assertEquals(1, mismatches.stream().filter(x -> x.MismatchType == MismatchType.REF_ONLY).count());
        assertEquals(1, mismatches.stream().filter(x -> x.MismatchType == MismatchType.NEW_ONLY).count());
        assertEquals(2, mismatches.stream().filter(x -> x.MismatchType == MismatchType.VALUE).count());
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
