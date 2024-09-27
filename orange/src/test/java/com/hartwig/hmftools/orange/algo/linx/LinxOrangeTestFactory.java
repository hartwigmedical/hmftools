package com.hartwig.hmftools.orange.algo.linx;

import com.hartwig.hmftools.common.linx.LinxTestFactory;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxBreakend;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxFusion;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxHomozygousDisruption;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxSvAnnotation;
import com.hartwig.hmftools.orange.conversion.LinxConversion;

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
                .from(LinxConversion.convert(LinxTestFactory.breakendBuilder().build()));
    }

    @NotNull
    public static ImmutableLinxHomozygousDisruption.Builder homozygousDisruptionBuilder()
    {
        return ImmutableLinxHomozygousDisruption.builder()
                .from(LinxConversion.convert(LinxTestFactory.homozygousDisruptionBuilder().build()));
    }
}
