package com.hartwig.hmftools.orange.algo.linx;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.linx.LinxDriver;
import com.hartwig.hmftools.common.linx.LinxFusion;
import com.hartwig.hmftools.common.linx.LinxGermlineDisruption;
import com.hartwig.hmftools.common.linx.LinxSvAnnotation;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface LinxData
{
    @NotNull
    List<LinxSvAnnotation> somaticSvAnnotations(); // used in interpreter to form final LinxBreakend info

    @NotNull
    List<LinxDriver> somaticDrivers();

    @NotNull
    List<LinxFusion> fusions();

    @NotNull
    List<LinxBreakend> somaticBreakends();

    @NotNull
    List<HomozygousDisruption> somaticHomozygousDisruptions();

    @Nullable
    List<LinxBreakend> germlineBreakends();

    @Nullable
    List<LinxGermlineDisruption> germlineDisruptions();

    @Nullable
    List<HomozygousDisruption> germlineHomozygousDisruptions();

    @Nullable
    List<LinxSvAnnotation> germlineSvAnnotations();

    @NotNull
    Set<Integer> fusionClusterIds();

    @NotNull
    Map<Integer, Integer> svIdToClusterId();

    @NotNull
    Map<Integer, Integer> clusterIdToLinkCount();

    @NotNull
    Map<Integer, Integer> clusterIdToExonCount();
}
