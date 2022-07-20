package com.hartwig.hmftools.cup;

import static com.hartwig.hmftools.cup.CuppaUtilsTest.TEST_CT_001;
import static com.hartwig.hmftools.cup.CuppaUtilsTest.TEST_CT_002;
import static com.hartwig.hmftools.cup.CuppaUtilsTest.TEST_SAMPLE_001;
import static com.hartwig.hmftools.cup.CuppaUtilsTest.TEST_SAMPLE_002;
import static com.hartwig.hmftools.cup.CuppaUtilsTest.TEST_SAMPLE_003;
import static com.hartwig.hmftools.common.cuppa.ClassifierType.GENOMIC_POSITION_COHORT;
import static com.hartwig.hmftools.common.cuppa.ClassifierType.SNV_96_PAIRWISE;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertTrue;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.cup.common.SampleData;
import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.common.SampleResult;
import com.hartwig.hmftools.cup.common.SampleSimilarity;
import com.hartwig.hmftools.cup.somatics.SomaticClassifier;

import org.junit.Test;

public class SomaticClassifierTest
{
    @Test
    public void testSnvCss()
    {
        SampleDataCache dataCache = new SampleDataCache();

        CuppaConfig config = new CuppaConfig();

        SomaticClassifier classifier = new SomaticClassifier(config, dataCache, null);

        CuppaUtilsTest.addRefSample(dataCache, TEST_SAMPLE_001, TEST_CT_001);
        CuppaUtilsTest.addRefSample(dataCache, TEST_SAMPLE_002, TEST_CT_002);

        final List<double[]> refSnvCounts = Lists.newArrayList();
        refSnvCounts.add(new double[] {20, 30, 40, 50});
        refSnvCounts.add(new double[] {20, 30, 42, 55});

        final List<double[]> refPosFreqCounts = Lists.newArrayList();
        refPosFreqCounts.add(new double[] {10, 12, 14, 16});
        refPosFreqCounts.add(new double[] {11, 12, 13, 15});

        final Map<String,double[]> cancerPosFreqCounts = Maps.newHashMap();
        cancerPosFreqCounts.put(TEST_CT_001, new double[] {20, 30, 40, 50});
        cancerPosFreqCounts.put(TEST_CT_002, new double[] {22, 30, 38, 52});

        classifier.addRefData(refSnvCounts, refPosFreqCounts, cancerPosFreqCounts);

        SampleData refSample1 = CuppaUtilsTest.addTestSample(dataCache, TEST_SAMPLE_001);
        SampleData refSample2 = CuppaUtilsTest.addTestSample(dataCache, TEST_SAMPLE_002);
        SampleData testSample = CuppaUtilsTest.addTestSample(dataCache, TEST_SAMPLE_003);

        final List<double[]> snvCounts = Lists.newArrayList();
        snvCounts.addAll(refSnvCounts);
        snvCounts.add(new double[] {1, 2, 3, 4});

        final List<double[]> posFreqCounts = Lists.newArrayList();
        posFreqCounts.addAll(refPosFreqCounts);
        posFreqCounts.add(new double[] {10, 12, 14, 16}); // exact match with ref 1, which will be penalised by subtraction

        List<String> sampleIds = Lists.newArrayList(refSample1.Id, refSample2.Id, testSample.Id);

        classifier.addSampleData(sampleIds, snvCounts, posFreqCounts);

        List<SampleResult> results = Lists.newArrayList();
        List<SampleSimilarity> similarities = Lists.newArrayList();

        classifier.processSample(testSample, results, similarities);

        SampleResult result = results.stream().filter(x -> x.DataType.equals(SNV_96_PAIRWISE.toString())).findFirst().orElse(null);
        assertTrue(result != null);
        assertEquals(0.35, result.CancerTypeValues.get(TEST_CT_001), 0.01);
        assertEquals(0.65, result.CancerTypeValues.get(TEST_CT_002), 0.01);

        result = results.stream().filter(x -> x.DataType.equals(GENOMIC_POSITION_COHORT.toString())).findFirst().orElse(null);
        assertTrue(result != null);
        assertEquals(0.48, result.CancerTypeValues.get(TEST_CT_001), 0.01);
        assertEquals(0.52, result.CancerTypeValues.get(TEST_CT_002), 0.01);

        // now test again with each of the ref samples
        results.clear();
        classifier.processSample(refSample1, results, similarities);

        // cannot match it's own cancer type
        result = results.stream().filter(x -> x.DataType.equals(SNV_96_PAIRWISE.toString())).findFirst().orElse(null);
        assertTrue(result != null);
        assertEquals(1.0, result.CancerTypeValues.get(TEST_CT_002), 0.01);

        // shifted way towards the CT 2 owing to subtraction
        result = results.stream().filter(x -> x.DataType.equals(GENOMIC_POSITION_COHORT.toString())).findFirst().orElse(null);
        assertTrue(result != null);
        assertEquals(0.01, result.CancerTypeValues.get(TEST_CT_001), 0.01);
        assertEquals(0.99, result.CancerTypeValues.get(TEST_CT_002), 0.01);

        results.clear();
        classifier.processSample(refSample2, results, similarities);

        result = results.stream().filter(x -> x.DataType.equals(SNV_96_PAIRWISE.toString())).findFirst().orElse(null);
        assertTrue(result != null);
        assertEquals(1.0, result.CancerTypeValues.get(TEST_CT_001), 0.01);

        result = results.stream().filter(x -> x.DataType.equals(GENOMIC_POSITION_COHORT.toString())).findFirst().orElse(null);
        assertTrue(result != null);
        assertEquals(0.99, result.CancerTypeValues.get(TEST_CT_001), 0.01);
        assertEquals(0.01, result.CancerTypeValues.get(TEST_CT_002), 0.01);

    }
}
