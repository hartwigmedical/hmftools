package com.hartwig.hmftools.orange.report.finding;

import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.finding.ImmutableSmallVariant;
import com.hartwig.hmftools.datamodel.driver.ReportedStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;
import com.hartwig.hmftools.orange.algo.purple.TestPurpleVariantFactory;
import com.hartwig.hmftools.orange.algo.util.PurpleDriverTestFactory;

import org.jetbrains.annotations.NotNull;

public final class TestVariantEntryFactory
{
    @NotNull
    public static ImmutableSmallVariant.Builder builder(String gene)
    {
        PurpleVariant variant = TestPurpleVariantFactory.builder()
                .gene(gene)
                .build();

        return ImmutableSmallVariant.builder()
                .findingKey("finding")
                .reportedStatus(ReportedStatus.REPORTED)
                .driverInterpretation(DriverInterpretation.HIGH)
                .purpleVariant(variant)
                .driver(PurpleDriverTestFactory.builder().build())
                .transcriptImpact(TestPurpleVariantFactory.impactBuilder().build())
                .isCanonical(true);
    }
}
