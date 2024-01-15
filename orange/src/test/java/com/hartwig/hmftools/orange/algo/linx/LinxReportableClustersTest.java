package com.hartwig.hmftools.orange.algo.linx;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.hartwig.hmftools.common.linx.DriverEventType;
import com.hartwig.hmftools.common.linx.ImmutableLinxData;
import com.hartwig.hmftools.common.linx.ImmutableLinxDriver;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.linx.LinxData;
import com.hartwig.hmftools.common.linx.LinxDriver;
import com.hartwig.hmftools.common.linx.LinxTestFactory;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class LinxReportableClustersTest
{
    @Test
    public void shouldFindAllReportableClusters()
    {
        LinxData linxData = linxDataBuilder().build();
        Set<Integer> clusters = LinxReportableClusters.findReportableClusters(linxData);
        assertEquals(new HashSet<>(Set.of(10, 50, 100)), clusters);
    }

    @Test
    public void shouldNotCountSimpleClusters()
    {
        // cluster 10 has single chain link and no exons -> should not be counted
        LinxData linxData = linxDataBuilder()
                .clusterIdToChainCount(Map.of(10, 1, 50, 2, 100, 2))
                .clusterIdToExonCount(Map.of(50, 1, 100, 1))
                .build();
        Set<Integer> clusters = LinxReportableClusters.findReportableClusters(linxData);
        assertEquals(new HashSet<>(Set.of(50, 100)), clusters);
    }

    @Test
    public void shouldCountMultiLinkClustersWithNoExon()
    {
        // cluster 10 has two chain links and no exons -> not simple so should be counted
        LinxData linxData2 = linxDataBuilder()
                .clusterIdToChainCount(Map.of(10, 2, 50, 2, 100, 2))
                .clusterIdToExonCount(Map.of(50, 1, 100, 1))
                .build();
        Set<Integer> clusters2 = LinxReportableClusters.findReportableClusters(linxData2);
        assertEquals(new HashSet<>(Set.of(10, 50, 100)), clusters2);
    }

    @Test
    public void shouldCountSingleLinkClusterWithExon()
    {
        // cluster 10 has single chain link and one exon -> not simple so should be counted
        LinxData linxData3 = linxDataBuilder()
                .clusterIdToChainCount(Map.of(10, 1, 50, 2, 100, 2))
                .clusterIdToExonCount(Map.of(10, 1, 50, 1, 100, 1))
                .build();
        Set<Integer> clusters3 = LinxReportableClusters.findReportableClusters(linxData3);
        assertEquals(new HashSet<>(Set.of(10, 50, 100)), clusters3);
    }

    @NotNull
    private static ImmutableLinxData.Builder linxDataBuilder()
    {
        List<LinxBreakend> breakends = new ArrayList<>();
        breakends.add(LinxTestFactory.breakendBuilder().svId(15).reportedDisruption(true).build());

        List<LinxDriver> drivers = new ArrayList<>();
        drivers.add(ImmutableLinxDriver.builder().clusterId(10).gene("Gene A").eventType(DriverEventType.DEL).build());

        List<Integer> fusions = List.of(50);
        Map<Integer, Integer> svToCluster = Map.of(15, 100);
        Map<Integer, Integer> clusterIdToChainCount = Map.of(10, 2, 50, 2, 100, 2);
        Map<Integer, Integer> clusterIdToExonCount = Map.of(10, 1, 50, 1, 100, 1);

        return ImmutableLinxData.builder()
                .reportableSomaticBreakends(breakends)
                .somaticDrivers(drivers)
                .fusionClusterIds(fusions)
                .svIdToClusterId(svToCluster)
                .clusterIdToChainCount(clusterIdToChainCount)
                .clusterIdToExonCount(clusterIdToExonCount);
    }
}
