package com.hartwig.hmftools.cup;

import static com.hartwig.hmftools.cup.CuppaConfig.formSamplePath;
import static com.hartwig.hmftools.common.cuppa.CategoryType.FEATURE;
import static com.hartwig.hmftools.cup.common.CupCalcs.adjustRefCounts;
import static com.hartwig.hmftools.cup.common.CupCalcs.calcCombinedFeatureResult;
import static com.hartwig.hmftools.cup.common.CupCalcs.convertToPercentages;
import static com.hartwig.hmftools.cup.common.CupConstants.CANCER_TYPE_UNKNOWN;
import static com.hartwig.hmftools.common.cuppa.ResultType.LIKELIHOOD;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertTrue;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.cuppa.ClassifierType;
import com.hartwig.hmftools.cup.common.SampleData;
import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.common.SampleResult;

import org.junit.Test;

public class CuppaUtilsTest
{
    public static final String TEST_SAMPLE_001 = "SAMPLE_001";
    public static final String TEST_SAMPLE_002 = "SAMPLE_002";
    public static final String TEST_SAMPLE_003 = "SAMPLE_003";
    public static final String TEST_SAMPLE_004 = "SAMPLE_004";

    public static final String TEST_CT_001 = "CT_001";
    public static final String TEST_CT_002 = "CT_002";
    public static final String TEST_CT_003 = "CT_003";
    public static final String TEST_CT_004 = "CT_004";

    public static SampleData addRefSample(final SampleDataCache dataCache, final String sampleId, final String cancerType)
    {
        SampleData sample = new SampleData(sampleId, cancerType, "");

        if(dataCache.findRefSampleData(sampleId) != null)
            return null;

        dataCache.addRefSample(sample);

        return sample;
    }

    public static SampleData addTestSample(final SampleDataCache dataCache, final String sampleId)
    {
        String cancerType = CANCER_TYPE_UNKNOWN;

        SampleData sample = new SampleData(sampleId, cancerType, "");

        if(dataCache.findSampleData(sampleId) != null)
            return null;

        dataCache.addTestSample(sample);

        return sample;
    }

    @Test
    public void testAdjustRefCalcs()
    {
        final double[] refCounts = new double[] {10, 20, 30, 40};
        final double[] sampleCounts = new double[] {2, 4, 6, 8};
        double[] adjustedCounts = adjustRefCounts(refCounts, sampleCounts, 1);
        assertEquals(8, adjustedCounts[0], 0.01);
        assertEquals(16, adjustedCounts[1], 0.01);
        assertEquals(24, adjustedCounts[2], 0.01);
        assertEquals(32, adjustedCounts[3], 0.01);

        adjustedCounts = adjustRefCounts(refCounts, sampleCounts, 0.5);
        assertEquals(9, adjustedCounts[0], 0.01);
        assertEquals(18, adjustedCounts[1], 0.01);
        assertEquals(27, adjustedCounts[2], 0.01);
        assertEquals(36, adjustedCounts[3], 0.01);
    }

    @Test
    public void testConfig()
    {
        String sampleWildcardPath = "/data/samples/*/";
        String sampleId = "TEST001";
        String samplePath = formSamplePath(sampleWildcardPath, sampleId);
        assertTrue(samplePath.equals("/data/samples/TEST001/"));
        assertTrue(sampleWildcardPath.equals("/data/samples/*/"));
    }

    @Test
    public void testConvertToPercentages()
    {
        // convertToPercentages
        Map<String,Double> dataMap = Maps.newHashMap();
        dataMap.put(TEST_CT_001, 0.1);
        dataMap.put(TEST_CT_002, 0.1);
        dataMap.put(TEST_CT_003, 0.1);
        dataMap.put(TEST_CT_004, 0.1);
        convertToPercentages(dataMap);

        assertEquals(0.25, dataMap.get(TEST_CT_001), 0.01);
        assertEquals(0.25, dataMap.get(TEST_CT_002), 0.01);
        assertEquals(0.25, dataMap.get(TEST_CT_003), 0.01);
        assertEquals(0.25, dataMap.get(TEST_CT_004), 0.01);

        dataMap = Maps.newHashMap();
        dataMap.put(TEST_CT_001, 1.0);
        dataMap.put(TEST_CT_002, 2.0);
        dataMap.put(TEST_CT_003, 3.0);
        dataMap.put(TEST_CT_004, 4.0);
        convertToPercentages(dataMap);

        assertEquals(0.1, dataMap.get(TEST_CT_001), 0.01);
        assertEquals(0.2, dataMap.get(TEST_CT_002), 0.01);
        assertEquals(0.3, dataMap.get(TEST_CT_003), 0.01);
        assertEquals(0.4, dataMap.get(TEST_CT_004), 0.01);
    }

    @Test
    public void testCombinedFeatureResult()
    {
        SampleData sample = new SampleData(TEST_SAMPLE_001, TEST_CT_001, "");

        List<SampleResult> results = Lists.newArrayList();

        Map<String,Double> cancerProbTotals = Maps.newHashMap();
        cancerProbTotals.put(TEST_CT_001, 0.5);
        cancerProbTotals.put(TEST_CT_002, 0.25);
        results.add(new SampleResult(TEST_SAMPLE_001, FEATURE, LIKELIHOOD, ClassifierType.FEATURE.toString(), "FEAT_01", cancerProbTotals));

        cancerProbTotals = Maps.newHashMap();
        cancerProbTotals.put(TEST_CT_001, 0.5);
        cancerProbTotals.put(TEST_CT_002, 0.5);

        results.add(new SampleResult(TEST_SAMPLE_001, FEATURE, LIKELIHOOD, ClassifierType.FEATURE.toString(), "FEAT_02", cancerProbTotals));

        SampleResult combResult = calcCombinedFeatureResult(sample, results, false);

        assertEquals(FEATURE, combResult.Category);
        assertEquals(0.63, combResult.CancerTypeValues.get(TEST_CT_001), 0.01);
        assertEquals(0.37, combResult.CancerTypeValues.get(TEST_CT_002), 0.01);
    }

}
