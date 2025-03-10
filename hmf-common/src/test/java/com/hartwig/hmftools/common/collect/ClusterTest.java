package com.hartwig.hmftools.common.collect;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.collect.Cluster.clusterCount;

import static org.junit.Assert.assertEquals;

import java.util.Collections;
import java.util.List;
import java.util.function.BiFunction;

import com.google.common.collect.Lists;

import org.junit.Test;

public class ClusterTest
{
    private static final BiFunction<Integer, Integer, Double> DISTANCE_FN = (a, b) -> (double) abs(a - b);

    @Test
    public void testClusterCountNoData()
    {
        List<Integer> data = Collections.emptyList();
        int actualClusterCount = clusterCount(data, DISTANCE_FN, 1);
        int expectedClustedCount = 0;

        assertEquals(expectedClustedCount, actualClusterCount);
    }

    @Test
    public void testClusterCountSingleDataPoint()
    {
        List<Integer> data = Lists.newArrayList(0);
        int actualClusterCount = clusterCount(data, DISTANCE_FN, 1);
        int expectedClustedCount = 1;

        assertEquals(expectedClustedCount, actualClusterCount);
    }

    @Test
    public void testClusterCountTwoDataPointsSingleCluster()
    {
        List<Integer> data = Lists.newArrayList(0, 1);
        int actualClusterCount = clusterCount(data, DISTANCE_FN, 1);
        int expectedClustedCount = 1;

        assertEquals(expectedClustedCount, actualClusterCount);
    }

    @Test
    public void testClusterCountTwoDataPointsTwoClusters()
    {
        List<Integer> data = Lists.newArrayList(0, 2);
        int actualClusterCount = clusterCount(data, DISTANCE_FN, 1);
        int expectedClustedCount = 2;

        assertEquals(expectedClustedCount, actualClusterCount);
    }

    @Test
    public void testClusterCountThreeDataPointsTwoClusters()
    {
        List<Integer> data = Lists.newArrayList(0, 1, 2);
        int actualClusterCount = clusterCount(data, DISTANCE_FN, 1);
        int expectedClustedCount = 2;

        assertEquals(expectedClustedCount, actualClusterCount);
    }

    @Test
    public void testClusterCountTwoClustersLargeSeparation()
    {
        List<Integer> data = Lists.newArrayList(0, 1, 2, 10, 11);
        int actualClusterCount = clusterCount(data, DISTANCE_FN, 1);
        int expectedClustedCount = 3;

        assertEquals(expectedClustedCount, actualClusterCount);
    }
}
