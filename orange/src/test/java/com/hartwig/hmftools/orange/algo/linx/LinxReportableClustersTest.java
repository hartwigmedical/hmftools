package com.hartwig.hmftools.orange.algo.linx;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.linx.DriverEventType;
import com.hartwig.hmftools.common.linx.ImmutableLinxData;
import com.hartwig.hmftools.common.linx.ImmutableLinxDriver;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.linx.LinxData;
import com.hartwig.hmftools.common.linx.LinxDriver;
import com.hartwig.hmftools.common.linx.LinxTestFactory;

import org.junit.Test;

public class LinxReportableClustersTest
{
    private static final String LINX_TEST_DATA_DIR = Resources.getResource("linx_reportable_clusters").getPath();

    @Test
    public void shouldFindAllReportableClusters() throws IOException
    {
        LinxData linxData = testLinxData();
        Set<Integer> clusters = LinxReportableClusters.findReportableClusters(linxData, LINX_TEST_DATA_DIR, "tumor_sample");
        assertEquals(new HashSet<>(Set.of(10, 50, 100)), clusters);
    }

    private static LinxData testLinxData()
    {
        List<LinxBreakend> breakends = new ArrayList<>();
        breakends.add(LinxTestFactory.breakendBuilder().svId(15).reportedDisruption(true).build());

        List<LinxDriver> drivers = new ArrayList<>();
        drivers.add(ImmutableLinxDriver.builder().clusterId(10).gene("Gene A").eventType(DriverEventType.DEL).build());

        LinxData linxData = ImmutableLinxData.builder().reportableSomaticBreakends(breakends).somaticDrivers(drivers).build();
        return linxData;
    }
}
