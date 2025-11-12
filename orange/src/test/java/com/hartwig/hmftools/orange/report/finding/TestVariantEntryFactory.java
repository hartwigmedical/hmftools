package com.hartwig.hmftools.orange.report.finding;

import com.hartwig.hmftools.datamodel.finding.DriverInterpretation;
import com.hartwig.hmftools.datamodel.finding.ImmutableSmallVariant;
import com.hartwig.hmftools.datamodel.purple.PurpleDriver;
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
                .isReportable(true)
                .isCandidate(false)
                .driverInterpretation(DriverInterpretation.HIGH)
                .purpleVariant(variant)
                .driver(PurpleDriverTestFactory.builder().build())
                .transcriptImpact(TestPurpleVariantFactory.impactBuilder().build())
                .isCanonical(true);
    }
}
