package com.hartwig.hmftools.cup;

import static com.hartwig.hmftools.cup.CuppaUtilsTest.TEST_CT_001;
import static com.hartwig.hmftools.cup.CuppaUtilsTest.TEST_CT_002;
import static com.hartwig.hmftools.cup.CuppaUtilsTest.TEST_SAMPLE_001;
import static com.hartwig.hmftools.cup.CuppaUtilsTest.TEST_SAMPLE_002;
import static com.hartwig.hmftools.cup.CuppaUtilsTest.TEST_SAMPLE_003;
import static com.hartwig.hmftools.cup.common.CategoryType.CLASSIFIER;
import static com.hartwig.hmftools.cup.common.ClassifierType.GENOMIC_POSITION_SIMILARITY;
import static com.hartwig.hmftools.cup.common.ClassifierType.SNV_96_PAIRWISE_SIMILARITY;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertTrue;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.cup.CuppaConfig;
import com.hartwig.hmftools.cup.common.SampleData;
import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.common.SampleResult;
import com.hartwig.hmftools.cup.common.SampleSimilarity;
import com.hartwig.hmftools.cup.somatics.SomaticClassifier;

import org.apache.commons.compress.utils.Lists;
import org.junit.Test;

public class SomaticClassifierTest
{

    @Test
    public void testSnvCss()
    {
        SampleDataCache dataCache = new SampleDataCache();

        CuppaConfig config = new CuppaConfig();

        SomaticClassifier classifier = new SomaticClassifier(config, dataCache);

        SampleData refSample1 = CuppaUtilsTest.addRefSample(dataCache, TEST_SAMPLE_001, TEST_CT_001);
        SampleData refSample2 = CuppaUtilsTest.addRefSample(dataCache, TEST_SAMPLE_002, TEST_CT_002);

        final List<double[]> refSnvCounts = Lists.newArrayList();
        refSnvCounts.add(new double[] {10, 20, 30, 40});
        refSnvCounts.add(new double[] {50, 50, 50, 50});

        final List<double[]> refPosFreqCounts = Lists.newArrayList();
        refPosFreqCounts.add(new double[] {1, 2, 3, 4});
        refPosFreqCounts.add(new double[] {5, 5, 5, 5});

        final Map<String,double[]> cancerPosFreqCounts = Maps.newHashMap();
        cancerPosFreqCounts.put(TEST_CT_001, new double[] {10, 20, 30, 40});
        cancerPosFreqCounts.put(TEST_CT_002, new double[] {40, 30, 20, 10});

        classifier.addRefData(refSnvCounts, refPosFreqCounts, cancerPosFreqCounts);

        SampleData testSample = CuppaUtilsTest.addTestSample(dataCache, TEST_SAMPLE_003);

        final List<double[]> snvCounts = Lists.newArrayList();
        snvCounts.add(new double[] {1, 2, 3, 4});

        final List<double[]> posFreqCounts = Lists.newArrayList();
        posFreqCounts.add(new double[] {5, 5, 5, 5});

        classifier.addSampleData(snvCounts, posFreqCounts);

        List<SampleResult> results = Lists.newArrayList();
        List<SampleSimilarity> similarities = Lists.newArrayList();

        classifier.processSample(testSample, results, similarities);

        SampleResult result = results.stream().filter(x -> x.DataType.equals(SNV_96_PAIRWISE_SIMILARITY.toString())).findFirst().orElse(null);
        assertTrue(result != null);
        assertEquals(1.0, result.CancerTypeValues.get(TEST_CT_001), 0.01);
        assertEquals(0.0, result.CancerTypeValues.get(TEST_CT_002), 0.01);

        result = results.stream().filter(x -> x.DataType.equals(GENOMIC_POSITION_SIMILARITY.toString())).findFirst().orElse(null);
        assertTrue(result != null);
        assertEquals(0.5, result.CancerTypeValues.get(TEST_CT_001), 0.01);
        assertEquals(0.5, result.CancerTypeValues.get(TEST_CT_002), 0.01);
    }
}
