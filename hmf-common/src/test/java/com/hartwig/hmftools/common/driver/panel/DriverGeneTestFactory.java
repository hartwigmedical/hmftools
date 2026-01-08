package com.hartwig.hmftools.common.driver.panel;

import static com.hartwig.hmftools.common.driver.panel.DriverGeneGermlineReporting.NONE;
import static com.hartwig.hmftools.common.purple.PurpleCommon.DEFAULT_DRIVER_AMPLIFICATION_PLOIDY_RATIO;
import static com.hartwig.hmftools.common.purple.PurpleCommon.DEFAULT_DRIVER_HET_DELETION_THRESHOLD;

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
                .amplificationRatio(DEFAULT_DRIVER_AMPLIFICATION_PLOIDY_RATIO)
                .reportHetDeletion(false)
                .reportLoh(false)
                .hetDeletionThreshold(DEFAULT_DRIVER_HET_DELETION_THRESHOLD)
                .reportSomaticHotspot(false)
                .reportGermlineVariant(NONE)
                .reportGermlineHotspot(NONE)
                .likelihoodType(DriverCategory.ONCO)
                .reportGermlineDisruption(NONE)
                .reportGermlineDeletion(NONE)
                .reportPGX(false);
    }
}
