package com.hartwig.hmftools.orange.algo.linx;

import static java.util.stream.Collectors.toList;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.linx.LinxData;

import org.jetbrains.annotations.NotNull;

// TODO this is redundant with linx.visualizer.SampleData.findReportableClusters(),
//  remove when available via shared code.
public final class LinxReportableClusters
{
    @NotNull
    public static Set<Integer> findReportableClusters(@NotNull LinxData linxData)
            throws IOException
    {
        Set<Integer> clusterIds = Sets.newHashSet();

        clusterIds.addAll(linxData.fusionClusterIds());

        List<LinxBreakend> breakends = linxData.reportableSomaticBreakends();
        List<Integer> svIds = breakends.stream().map(LinxBreakend::svId).collect(toList());

        Map<Integer, Integer> svToCluster = linxData.SvIdToClusterId();
        svIds.stream().filter(svToCluster::containsKey).map(svToCluster::get).forEach(clusterIds::add);

        linxData.somaticDrivers().stream().filter(x -> x.clusterId() >= 0).forEach(x -> clusterIds.add(x.clusterId()));
        return clusterIds;
    }
}
