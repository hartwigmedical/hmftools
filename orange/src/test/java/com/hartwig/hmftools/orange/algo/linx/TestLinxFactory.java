package com.hartwig.hmftools.orange.algo.linx;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.linx.ImmutableLinxBreakend;
import com.hartwig.hmftools.common.linx.ImmutableLinxFusion;
import com.hartwig.hmftools.common.linx.ImmutableLinxSvAnnotation;
import com.hartwig.hmftools.common.linx.LinxTestFactory;

public final class TestLinxFactory
{
    public static ImmutableLinxFusion.Builder fusionBuilder()
    {
        return LinxTestFactory.fusionBuilder();
    }

    public static ImmutableLinxBreakend.Builder breakendBuilder()
    {
        return LinxTestFactory.breakendBuilder();
    }

    public static ImmutableLinxSvAnnotation.Builder svAnnotationBuilder()
    {
        return LinxTestFactory.svAnnotationBuilder();
    }

    public static ImmutableLinxData.Builder linxDataBuilder()
    {
        return ImmutableLinxData.builder()
                .somaticDrivers(Lists.newArrayList())
                .somaticBreakends(Lists.newArrayList())
                .somaticSvAnnotations(Lists.newArrayList())
                .somaticHomozygousDisruptions(Lists.newArrayList())
                .somaticDriverData(Lists.newArrayList())
                .somaticBreakendClusterIds(Maps.newHashMap())
                .fusions(Lists.newArrayList())
                .fusionClusterIds(Maps.newHashMap())
                .germlineDisruptions(Lists.newArrayList())
                .germlineBreakends(Lists.newArrayList())
                .germlineDrivers(Lists.newArrayList())
                .germlineSvAnnotations(Lists.newArrayList())
                .reportableEventPlots(Lists.newArrayList());
    }
}
