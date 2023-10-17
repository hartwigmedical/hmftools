package com.hartwig.hmftools.orange.algo.linx;

import static java.util.stream.Collectors.toList;

import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.linx.LinxData;

import org.jetbrains.annotations.NotNull;

// TODO this is redundant with linx.visualizer.SampleData.findReportableClusters(),
//  remove when available via shared code.
public final class LinxReportableClusters
{
    @NotNull
    public static Set<Integer> findReportableClusters(@NotNull LinxData linx)
    {
        Set<Integer> clusterIds = new HashSet<>();

        clusterIds.addAll(linx.fusionClusterIds());

        List<LinxBreakend> breakends = linx.reportableSomaticBreakends();
        List<Integer> svIds = breakends.stream().map(LinxBreakend::svId).collect(toList());

        Map<Integer, Integer> svToCluster = linx.svIdToClusterId();
        svIds.stream().filter(svToCluster::containsKey).map(svToCluster::get).forEach(clusterIds::add);

        linx.somaticDrivers().stream().filter(x -> x.clusterId() >= 0).forEach(x -> clusterIds.add(x.clusterId()));
        return clusterIds;
    }
}
