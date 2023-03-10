package com.hartwig.hmftools.orange.algo.linx;

import com.hartwig.hmftools.common.linx.LinxTestFactory;
import com.hartwig.hmftools.datamodel.linx.*;
import org.jetbrains.annotations.NotNull;

public class LinxOrangeTestFactory {

    @NotNull
    public static ImmutableLinxSvAnnotation.Builder svAnnotationBuilder() {
        return ImmutableLinxSvAnnotation.builder()
                .from(LinxInterpreter.convert(LinxTestFactory.svAnnotationBuilder().build()));
    }

    @NotNull
    public static ImmutableLinxFusion.Builder fusionBuilder() {
        return ImmutableLinxFusion.builder()
                .from(LinxInterpreter.convert(LinxTestFactory.fusionBuilder().build()));
    }

    @NotNull
    public static ImmutableLinxBreakend.Builder breakendBuilder() {
        return ImmutableLinxBreakend.builder()
                .from(LinxInterpreter.convert(LinxTestFactory.breakendBuilder().build()));
    }

    @NotNull
    public static ImmutableHomozygousDisruption.Builder homozygousDisruptionBuilder() {
        return ImmutableHomozygousDisruption.builder()
                .from(LinxInterpreter.convert(LinxTestFactory.homozygousDisruptionBuilder().build()));
    }
}
