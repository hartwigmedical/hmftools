package com.hartwig.hmftools.cobalt.norm;

import static com.hartwig.hmftools.cobalt.norm.DataLoader.addTargetRegions;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.bed.ImmutableNamedBed;
import com.hartwig.hmftools.common.genome.bed.NamedBed;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

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
        List<NamedBed> namedBedRecords = Lists.newArrayList();

        namedBedRecords.add(ImmutableNamedBed.builder().name("REGION1").chromosome(CHR_1).start(1100).end(2900).build());
        namedBedRecords.add(ImmutableNamedBed.builder().name("REGION2").chromosome(CHR_1).start(2950).end(4000).build());

        namedBedRecords.add(ImmutableNamedBed.builder().name("REGION3").chromosome(CHR_1).start(5000).end(6001).build());

        namedBedRecords.add(ImmutableNamedBed.builder().name("REGION3").chromosome(CHR_2).start(9000).end(10001).build());

        namedBedRecords.add(ImmutableNamedBed.builder().name("REGION3").chromosome(CHR_2).start(12999).end(14500).build());
        Map<String,List<RegionData>> chrRegionData = Maps.newHashMap();

        addTargetRegions(namedBedRecords, chrRegionData);

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
