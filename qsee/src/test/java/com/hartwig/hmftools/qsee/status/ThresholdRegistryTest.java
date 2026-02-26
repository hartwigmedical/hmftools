package com.hartwig.hmftools.qsee.status;

import static com.hartwig.hmftools.qsee.common.QseeConstants.QC_LOGGER;

import static org.junit.Assert.assertEquals;

import com.google.common.io.Resources;
import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.FeatureType;
import com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
import org.junit.Test;

public class ThresholdRegistryTest
{
    private static final String TEST_THRESHOLD_OVERRIDES_FILE = Resources.getResource("threshold_overrides.tsv").getPath();

    public ThresholdRegistryTest(){ Configurator.setLevel(QC_LOGGER.getName(), Level.DEBUG); }

    @Test
    public void canOverrideThresholdsFromFile()
    {
        ThresholdRegistry thresholds = ThresholdOverridesFile.read(TEST_THRESHOLD_OVERRIDES_FILE);

        ThresholdKey key;
        QcThreshold expectedThreshold;

        key = new ThresholdKey(SampleType.TUMOR, FeatureType.SUMMARY_TABLE, SummaryTableFeature.MAPPED_PROPORTION.name(), QcStatusType.FAIL);
        expectedThreshold = new QcThreshold(key, ComparisonOperator.LESS_THAN, Double.NEGATIVE_INFINITY);
        assertEquals(thresholds.getThreshold(key), expectedThreshold);

        key = new ThresholdKey(SampleType.TUMOR, FeatureType.SUMMARY_TABLE, SummaryTableFeature.LOW_MAP_QUAL.name(), QcStatusType.WARN);
        expectedThreshold = new QcThreshold(key, ComparisonOperator.GREATER_THAN, Double.NaN);
        assertEquals(thresholds.getThreshold(key), expectedThreshold);
    }
}
