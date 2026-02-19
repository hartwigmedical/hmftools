package com.hartwig.hmftools.orange.report.datamodel;

import com.hartwig.hmftools.datamodel.purple.HotspotType;
import com.hartwig.hmftools.datamodel.purple.PurpleGenotypeStatus;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class TestVariantEntryFactory
{
    @NotNull
    public static ImmutableVariantEntry.Builder builder()
    {
        return ImmutableVariantEntry.builder()
                .gene(Strings.EMPTY)
                .isCanonical(false)
                .impact(Strings.EMPTY)
                .affectedCodon(1)
                .variantCopyNumber(0D)
                .totalCopyNumber(0D)
                .minorAlleleCopyNumber(0D)
                .biallelic(false)
                .biallelicProbability(0.1)
                .hotspot(HotspotType.NON_HOTSPOT)
                .driverLikelihood(1.0)
                .clonalLikelihood(0D)
                .somaticLikelihood("")
                .vaf(0.5)
                .depth(10)
                .rnaDepth(null)
                .genotypeStatus(PurpleGenotypeStatus.UNKNOWN);
    }
}
