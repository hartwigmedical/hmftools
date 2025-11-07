package com.hartwig.hmftools.common.driver;

import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.ReportedStatus;

public final class DriverCatalogFactory
{
    private DriverCatalogFactory() {}

    @Deprecated
    public static DriverCatalog createCopyNumberDriver(
            DriverCategory category, DriverType driver, final LikelihoodMethod likelihoodMethod, final boolean biallelic,
            final double likelihood, final GeneCopyNumber geneCopyNumber)
    {
        return ImmutableDriverCatalog.builder()
                .chromosome(geneCopyNumber.chromosome())
                .chromosomeBand(geneCopyNumber.ChromosomeBand)
                .gene(geneCopyNumber.geneName())
                .transcript(geneCopyNumber.TransName)
                .isCanonical(geneCopyNumber.IsCanonical)
                .missense(0)
                .nonsense(0)
                .inframe(0)
                .frameshift(0)
                .splice(0)
                .driverLikelihood(likelihood)
                .driver(driver)
                .likelihoodMethod(likelihoodMethod)
                .category(category)
                .reportedStatus(ReportedStatus.REPORTED)
                .biallelic(biallelic)
                .minCopyNumber(geneCopyNumber.minCopyNumber())
                .maxCopyNumber(geneCopyNumber.maxCopyNumber())
                .reportedStatus(ReportedStatus.NONE)
                .build();
    }

}
