package com.hartwig.hmftools.cobalt.norm;

import static com.hartwig.hmftools.cobalt.norm.DataLoader.addTargetRegions;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.bed.ImmutableNamedBed;
import com.hartwig.hmftools.common.genome.bed.NamedBed;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import static com.hartwig.hmftools.cobalt.norm.Normaliser.calcSampleAdjustedRatios;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class NormalisationFileTest
{
    @Test
    public void testCreateTargetRegions()
    {
        List<ChrBaseRegion> bedRegions = Lists.newArrayList();

        bedRegions.add(new ChrBaseRegion(CHR_1, 1100, 2900));
        bedRegions.add(new ChrBaseRegion(CHR_1, 2950, 4000));
        bedRegions.add(new ChrBaseRegion(CHR_1, 5000, 6001));
        bedRegions.add(new ChrBaseRegion(CHR_2, 9000, 10001));
        bedRegions.add(new ChrBaseRegion(CHR_2, 12999, 14500));

        Map<String,List<RegionData>> chrRegionData = Maps.newHashMap();

        addTargetRegions(bedRegions, chrRegionData);

        List<RegionData> allRegions = Lists.newArrayList();
        chrRegionData.values().forEach(x -> allRegions.addAll(x));

        assertEquals(10, allRegions.size());
        int index = 0;
        assertEquals(1001, allRegions.get(index++).Position);
        assertEquals(2001, allRegions.get(index++).Position);

        assertEquals(3001, allRegions.get(index++).Position);

        assertEquals(5001, allRegions.get(index++).Position);
        assertEquals(6001, allRegions.get(index++).Position);

        assertEquals(9001, allRegions.get(index++).Position);
        assertEquals(10001, allRegions.get(index++).Position);
        assertEquals(12001, allRegions.get(index++).Position);
        assertEquals(13001, allRegions.get(index++).Position);
        assertEquals(14001, allRegions.get(index++).Position);
    }

    @Test
    public void testSampleNormalisation()
    {
        Map<String,List<RegionData>> chrRegionData = Maps.newHashMap();

        List<RegionData> regions = Lists.newArrayList();
        regions.add(createRegion(1001, 15, 1.0, 10, 0.5)); // filtered
        regions.add(createRegion(2001, 20, 1.0, 20, 0.5));
        regions.add(createRegion(3001, 30, 1.0, 30, 0.5));
        regions.add(createRegion(3001, 30, 1.0, -1, 0.5)); // filtered
        regions.add(createRegion(4001, 40, 1.0, 40, 0.5));
        regions.add(createRegion(4001, 50, 0.5, 50, 0.5)); // filtered
        regions.add(createRegion(4001, 60, 1.0, 60, 0.5));
        regions.add(createRegion(5001, 70, 1.0, 70, 0.5)); // filtered

        chrRegionData.put(CHR_1, regions);
        chrRegionData.put(HumanChromosome._X.toString(), regions);
        chrRegionData.put(HumanChromosome._Y.toString(), regions);

        NormCalcData normCalcData = calcSampleAdjustedRatios(0, chrRegionData);
        assertEquals(4, normCalcData.SampleFilteredRegionCount);
        assertEquals(37.5, normCalcData.SampleMeanReadDepth, 0.1);
        assertEquals(35, normCalcData.SampleMedianReadDepth, 0.1);
        assertEquals(4, normCalcData.GcBucketMedians.size());
    }

    @SuppressWarnings("SameParameterValue")
    private static RegionData createRegion(int position, int gcBucket, double mappability, int readCount, double gcRatioPanel)
    {
        RegionData regionData = new RegionData(position);
        regionData.setGcProfile(gcBucket, mappability);
        regionData.addSampleRegionData(new SampleRegionData(readCount, gcRatioPanel, 1));
        return regionData;
    }

}
