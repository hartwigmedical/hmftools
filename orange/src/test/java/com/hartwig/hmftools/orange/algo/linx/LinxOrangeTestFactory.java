package com.hartwig.hmftools.orange.algo.linx;

import com.hartwig.hmftools.common.linx.LinxTestFactory;
import com.hartwig.hmftools.datamodel.gene.TranscriptCodingType;
import com.hartwig.hmftools.datamodel.gene.TranscriptRegionType;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxBreakend;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxFusion;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxHomozygousDisruption;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxSvAnnotation;
import com.hartwig.hmftools.datamodel.linx.LinxBreakendType;
import com.hartwig.hmftools.orange.conversion.LinxConversion;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class LinxOrangeTestFactory
{
    @NotNull
    public static ImmutableLinxSvAnnotation.Builder svAnnotationBuilder()
    {
        return ImmutableLinxSvAnnotation.builder()
                .from(LinxConversion.convert(LinxTestFactory.svAnnotationBuilder().build()));
    }

    @NotNull
    public static ImmutableLinxFusion.Builder fusionBuilder()
    {
        return ImmutableLinxFusion.builder()
                .from(LinxConversion.convert(LinxTestFactory.fusionBuilder().build()));
    }

    @NotNull
    public static ImmutableLinxBreakend.Builder breakendBuilder()
    {
        return ImmutableLinxBreakend.builder()
                .id(0)
                .svId(0)
                .gene("GENE")
                .chromosome("1")
                .chromosomeBand("2")
                .transcript(Strings.EMPTY)
                .isCanonical(true)
                .geneOrientation(Strings.EMPTY)
                .disruptive(false)
                .reported(false)
                .undisruptedCopyNumber(0D)
                .type(LinxBreakendType.BND)
                .regionType(TranscriptRegionType.UNKNOWN)
                .codingType(TranscriptCodingType.UNKNOWN)
                .nextSpliceExonRank(0)
                .orientation(0)
                .exonUp(0)
                .exonDown(0)
                .junctionCopyNumber(0D);
    }
}
