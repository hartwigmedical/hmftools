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
        LinxData linxData = createTestLinxData();
        Set<Integer> clusters = LinxReportableClusters.findReportableClusters(linxData);
        assertEquals(new HashSet<>(Set.of(10, 50, 100)), clusters);
    }

    @NotNull
    private static LinxData createTestLinxData()
    {
        List<LinxBreakend> breakends = new ArrayList<>();
        breakends.add(LinxTestFactory.breakendBuilder().svId(15).reportedDisruption(true).build());

        List<LinxDriver> drivers = new ArrayList<>();
        drivers.add(ImmutableLinxDriver.builder().clusterId(10).gene("Gene A").eventType(DriverEventType.DEL).build());

        List<Integer> fusions = List.of(50);
        Map<Integer, Integer> svToCluster = Map.of(15, 100);

        return ImmutableLinxData.builder()
                .reportableSomaticBreakends(breakends)
                .somaticDrivers(drivers)
                .fusionClusterIds(fusions)
                .putAllSvIdToClusterId(svToCluster)
                .build();
    }
}
