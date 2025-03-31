package com.hartwig.hmftools.common.driver.panel;

import static com.hartwig.hmftools.common.driver.panel.DriverGeneGermlineReporting.NONE;

import com.hartwig.hmftools.common.driver.DriverCategory;

import org.apache.logging.log4j.util.Strings;

public final class DriverGeneTestFactory
{
    public static ImmutableDriverGene.Builder builder()
    {
        return ImmutableDriverGene.builder()
                .gene(Strings.EMPTY)
                .reportMissenseAndInframe(false)
                .reportNonsenseAndFrameshift(false)
                .reportSplice(false)
                .reportDeletion(false)
                .reportDisruption(false)
                .reportAmplification(false)
                .reportSomaticHotspot(false)
                .reportGermlineVariant(NONE)
                .reportGermlineHotspot(NONE)
                .likelihoodType(DriverCategory.ONCO)
                .reportGermlineDisruption(NONE)
                .reportGermlineDeletion(NONE)
                .reportPGX(false);
    }
}
