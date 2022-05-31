package com.hartwig.hmftools.cup;

import static com.hartwig.hmftools.cup.CuppaUtilsTest.TEST_CT_001;
import static com.hartwig.hmftools.cup.CuppaUtilsTest.TEST_CT_002;
import static com.hartwig.hmftools.cup.CuppaUtilsTest.TEST_CT_003;
import static com.hartwig.hmftools.cup.CuppaUtilsTest.TEST_SAMPLE_001;
import static com.hartwig.hmftools.cup.CuppaUtilsTest.TEST_SAMPLE_002;
import static com.hartwig.hmftools.cup.CuppaUtilsTest.TEST_SAMPLE_003;
import static com.hartwig.hmftools.cup.common.ClassifierType.FEATURE;
import static com.hartwig.hmftools.cup.feature.FeatureType.AMP;
import static com.hartwig.hmftools.cup.feature.FeatureType.DRIVER;
import static com.hartwig.hmftools.cup.feature.FeatureType.FUSION;
import static com.hartwig.hmftools.cup.feature.FeatureType.VIRUS;

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
import com.hartwig.hmftools.cup.feature.FeatureClassifier;
import com.hartwig.hmftools.cup.feature.FeaturePrevData;
import com.hartwig.hmftools.cup.feature.SampleFeatureData;

import org.junit.Test;

public class FeatureClassifierTest
{
    private static final String FUSION_001 = "FUSION_001";
    private static final String VIRUS_001 = "VIRUS_001";
    private static final String GENE_001 = "GENE_001";
    private static final String GENE_002 = "GENE_002";

    @Test
    public void testFeaturesNonRefSamples()
    {
        SampleDataCache dataCache = new SampleDataCache();

        CuppaConfig config = new CuppaConfig();

        FeatureClassifier classifier = new FeatureClassifier(config, dataCache, null);

        final Map<String,List<SampleFeatureData>> sampleFeatureMap = Maps.newHashMap();
        final Map<String,List<FeaturePrevData>> cancerFeaturePrevalence = Maps.newHashMap();

        List<FeaturePrevData> prevalences = Lists.newArrayList();
        prevalences.add(new FeaturePrevData(TEST_CT_001, FUSION_001, FUSION, 0.4));
        prevalences.add(new FeaturePrevData(TEST_CT_001, GENE_001, DRIVER, 0.1));
        cancerFeaturePrevalence.put(TEST_CT_001, prevalences);

        prevalences = Lists.newArrayList();
        prevalences.add(new FeaturePrevData(TEST_CT_002, FUSION_001, FUSION, 0.1));
        prevalences.add(new FeaturePrevData(TEST_CT_002, GENE_002, DRIVER, 0.6));
        cancerFeaturePrevalence.put(TEST_CT_002, prevalences);

        prevalences = Lists.newArrayList();
        prevalences.add(new FeaturePrevData(TEST_CT_003, GENE_002, DRIVER, 0.1));
        prevalences.add(new FeaturePrevData(TEST_CT_003, "VIRUS_001", VIRUS, 0.2));
        cancerFeaturePrevalence.put(TEST_CT_003, prevalences);

        SampleData testSample1 = CuppaUtilsTest.addTestSample(dataCache, TEST_SAMPLE_001);
        SampleData testSample2 = CuppaUtilsTest.addTestSample(dataCache, TEST_SAMPLE_002);
        SampleData testSample3 = CuppaUtilsTest.addTestSample(dataCache, TEST_SAMPLE_003);

        List<SampleFeatureData> sampleFeatures = Lists.newArrayList();
        sampleFeatures.add(new SampleFeatureData(testSample1.Id, FUSION_001, FUSION, 1));
        sampleFeatureMap.put(testSample1.Id, sampleFeatures);

        sampleFeatures = Lists.newArrayList();
        sampleFeatures.add(new SampleFeatureData(testSample2.Id, GENE_001, DRIVER, 1));
        sampleFeatureMap.put(testSample2.Id, sampleFeatures);

        sampleFeatures = Lists.newArrayList();
        sampleFeatures.add(new SampleFeatureData(testSample3.Id, VIRUS_001, VIRUS, 1));
        sampleFeatures.add(new SampleFeatureData(testSample3.Id, GENE_002, DRIVER, 1));
        sampleFeatureMap.put(testSample3.Id, sampleFeatures);

        classifier.addFeaturePrevalences(sampleFeatureMap, cancerFeaturePrevalence);

        List<SampleResult> results = Lists.newArrayList();
        List<SampleSimilarity> similarities = Lists.newArrayList();

        classifier.processSample(testSample1, results, similarities);

        SampleResult result = results.stream().filter(x -> x.DataType.equals(FEATURE.toString())).findFirst().orElse(null);
        assertTrue(result != null);
        assertEquals(0.79, result.CancerTypeValues.get(TEST_CT_001), 0.01);
        assertEquals(0.21, result.CancerTypeValues.get(TEST_CT_002), 0.01);
        assertEquals(0.01, result.CancerTypeValues.get(TEST_CT_003), 0.01);

        results.clear();
        classifier.processSample(testSample2, results, similarities);

        result = results.stream().filter(x -> x.DataType.equals(FEATURE.toString())).findFirst().orElse(null);
        assertTrue(result != null);
        assertEquals(0.67, result.CancerTypeValues.get(TEST_CT_001), 0.01);
        assertEquals(0.17, result.CancerTypeValues.get(TEST_CT_002), 0.01);
        assertEquals(0.17, result.CancerTypeValues.get(TEST_CT_003), 0.01);

        results.clear();
        classifier.processSample(testSample3, results, similarities);

        result = results.stream().filter(x -> x.DataType.equals(FEATURE.toString())).findFirst().orElse(null);
        assertTrue(result != null);
        assertEquals(0.0, result.CancerTypeValues.get(TEST_CT_001), 0.01);
        assertEquals(0.02, result.CancerTypeValues.get(TEST_CT_002), 0.01);
        assertEquals(0.16, result.CancerTypeValues.get(TEST_CT_003), 0.01);
    }

    @Test
    public void testDriverAmps()
    {
        SampleDataCache dataCache = new SampleDataCache();

        CuppaConfig config = new CuppaConfig();

        FeatureClassifier classifier = new FeatureClassifier(config, dataCache, null);

        final Map<String, List<SampleFeatureData>> sampleFeatureMap = Maps.newHashMap();
        final Map<String, List<FeaturePrevData>> cancerFeaturePrevalence = Maps.newHashMap();

        List<FeaturePrevData> prevalences = Lists.newArrayList();
        prevalences.add(new FeaturePrevData(TEST_CT_001, GENE_001, DRIVER, 0.2));
        prevalences.add(new FeaturePrevData(TEST_CT_001, GENE_001, AMP, 0.8));
        cancerFeaturePrevalence.put(TEST_CT_001, prevalences);

        prevalences = Lists.newArrayList();
        prevalences.add(new FeaturePrevData(TEST_CT_002, GENE_001, DRIVER, 0.8));
        prevalences.add(new FeaturePrevData(TEST_CT_002, GENE_001, AMP, 0.2));
        cancerFeaturePrevalence.put(TEST_CT_002, prevalences);

        SampleData testSample1 = CuppaUtilsTest.addTestSample(dataCache, TEST_SAMPLE_001);
        SampleData testSample2 = CuppaUtilsTest.addTestSample(dataCache, TEST_SAMPLE_002);

        List<SampleFeatureData> sampleFeatures = Lists.newArrayList();
        sampleFeatures.add(new SampleFeatureData(testSample1.Id, GENE_001, AMP, 1));
        sampleFeatureMap.put(testSample1.Id, sampleFeatures);

        sampleFeatures = Lists.newArrayList();
        sampleFeatures.add(new SampleFeatureData(testSample2.Id, GENE_001, DRIVER, 1));
        sampleFeatureMap.put(testSample2.Id, sampleFeatures);

        classifier.addFeaturePrevalences(sampleFeatureMap, cancerFeaturePrevalence);

        List<SampleResult> results = Lists.newArrayList();
        List<SampleSimilarity> similarities = Lists.newArrayList();

        classifier.processSample(testSample1, results, similarities);

        SampleResult result = results.stream().filter(x -> x.DataType.equals(FEATURE.toString())).findFirst().orElse(null);
        assertTrue(result != null);
        assertEquals(0.77, result.CancerTypeValues.get(TEST_CT_001), 0.01);
        assertEquals(0.23, result.CancerTypeValues.get(TEST_CT_002), 0.01);

        results.clear();
        classifier.processSample(testSample2, results, similarities);

        result = results.stream().filter(x -> x.DataType.equals(FEATURE.toString())).findFirst().orElse(null);
        assertTrue(result != null);
        assertEquals(0.23, result.CancerTypeValues.get(TEST_CT_001), 0.01);
        assertEquals(0.77, result.CancerTypeValues.get(TEST_CT_002), 0.01);
    }

    @Test
    public void testFeaturesRefSamples()
    {
        SampleDataCache dataCache = new SampleDataCache();

        CuppaConfig config = new CuppaConfig();

        FeatureClassifier classifier = new FeatureClassifier(config, dataCache, null);

        final Map<String,List<SampleFeatureData>> sampleFeatureMap = Maps.newHashMap();
        final Map<String,List<FeaturePrevData>> cancerFeaturePrevalence = Maps.newHashMap();

        List<FeaturePrevData> prevalences = Lists.newArrayList();
        prevalences.add(new FeaturePrevData(TEST_CT_001, GENE_001, DRIVER, 0.4));
        cancerFeaturePrevalence.put(TEST_CT_001, prevalences);

        prevalences = Lists.newArrayList();
        prevalences.add(new FeaturePrevData(TEST_CT_002, GENE_002, DRIVER, 0.3));
        cancerFeaturePrevalence.put(TEST_CT_002, prevalences);

        // feature incidences need to be supported by actual ref samples to allow for ref-sample adjustments
        int refSampleCount = 10;
        int refSampleIndex = 0;
        for(int i = 0; i < refSampleCount; ++i)
        {
            CuppaUtilsTest.addRefSample(dataCache, String.format("SAMPLE_%03d", refSampleIndex), TEST_CT_001);
            ++refSampleIndex;
        }

        for(int i = 0; i < refSampleCount; ++i)
        {
            CuppaUtilsTest.addRefSample(dataCache, String.format("SAMPLE_%03d", refSampleIndex), TEST_CT_002);
            ++refSampleIndex;
        }

        String refSampleId = dataCache.RefCancerSampleData.get(TEST_CT_001).get(0).Id;
        SampleData refSample1 = CuppaUtilsTest.addTestSample(dataCache, refSampleId);

        refSampleId = dataCache.RefCancerSampleData.get(TEST_CT_002).get(0).Id;
        SampleData refSample2 = CuppaUtilsTest.addTestSample(dataCache, refSampleId);

        List<SampleFeatureData> sampleFeatures = Lists.newArrayList();
        sampleFeatures.add(new SampleFeatureData(refSample1.Id, GENE_001, DRIVER, 1));
        sampleFeatureMap.put(refSample1.Id, sampleFeatures);

        sampleFeatures = Lists.newArrayList();
        sampleFeatures.add(new SampleFeatureData(refSample2.Id, GENE_002, DRIVER, 1));
        sampleFeatureMap.put(refSample2.Id, sampleFeatures);

        // also add a test sample with the same feature as the ref sample to highlight the discounting of ref prevalence
        SampleData testSample1 = CuppaUtilsTest.addTestSample(dataCache, "SAMPLE_TEST");

        sampleFeatures = Lists.newArrayList();
        sampleFeatures.add(new SampleFeatureData(testSample1.Id, GENE_001, DRIVER, 1));
        sampleFeatureMap.put(testSample1.Id, sampleFeatures);

        classifier.addFeaturePrevalences(sampleFeatureMap, cancerFeaturePrevalence);

        List<SampleResult> results = Lists.newArrayList();
        List<SampleSimilarity> similarities = Lists.newArrayList();

        classifier.processSample(refSample1, results, similarities);

        SampleResult result = results.stream().filter(x -> x.DataType.equals(FEATURE.toString())).findFirst().orElse(null);
        assertTrue(result != null);
        assertEquals(0.88, result.CancerTypeValues.get(TEST_CT_001), 0.01);
        assertEquals(0.10, result.CancerTypeValues.get(TEST_CT_002), 0.01);

        // compare non-ref sample with same feature - will have a higher result
        results.clear();
        classifier.processSample(testSample1, results, similarities);

        result = results.stream().filter(x -> x.DataType.equals(FEATURE.toString())).findFirst().orElse(null);
        assertTrue(result != null);
        assertEquals(0.9, result.CancerTypeValues.get(TEST_CT_001), 0.01);
        assertEquals(0.10, result.CancerTypeValues.get(TEST_CT_002), 0.01);

        results.clear();
        classifier.processSample(refSample2, results, similarities);

        result = results.stream().filter(x -> x.DataType.equals(FEATURE.toString())).findFirst().orElse(null);
        assertTrue(result != null);
        assertEquals(0.125, result.CancerTypeValues.get(TEST_CT_001), 0.01);
        assertEquals(0.84, result.CancerTypeValues.get(TEST_CT_002), 0.01);
    }

}
