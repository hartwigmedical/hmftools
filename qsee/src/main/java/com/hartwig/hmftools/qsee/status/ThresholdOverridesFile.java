package com.hartwig.hmftools.qsee.status;

import static com.hartwig.hmftools.qsee.common.QseeConstants.QC_LOGGER;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.COL_FEATURE_NAME;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.COL_FEATURE_TYPE;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.COL_SAMPLE_TYPE;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.utils.file.DelimFileReader;
import com.hartwig.hmftools.common.utils.file.DelimFileWriter;
import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.FeatureType;

public class ThresholdOverridesFile
{
    public static final String THRESHOLD_OVERRIDES_FILE_CFG = "threshold_overrides_file";
    public static final String THRESHOLD_OVERRIDES_FILE_CFG_DESC = "Path to threshold overrides file";

    private static final String COL_QC_STATUS_TYPE = "QcStatusType";
    private static final String COL_COMPARISON_OPERATOR = "ComparisonOperator";
    private static final String COL_THRESHOLD = "Threshold";
    private static final String COL_OVERRIDABLE = "Overridable";

    public static ThresholdRegistry read(String filename, boolean targetedMode)
    {
        try(DelimFileReader reader = new DelimFileReader(filename))
        {
            QC_LOGGER.info("Reading threshold overrides from: {}", filename);

            ThresholdRegistry thresholds = ThresholdRegistry.createEmptyModifiable().setDefaults(targetedMode);
            boolean hasUnsupportedThresholds = false;
            List<ThresholdKey> seenKeys = new ArrayList<>();

            for(DelimFileReader.Row row : reader)
            {
                SampleType sampleType = SampleType.valueOf(row.get(COL_SAMPLE_TYPE));
                FeatureType featureType = FeatureType.valueOf(row.get(COL_FEATURE_TYPE));
                String featureName = row.get(COL_FEATURE_NAME);
                QcStatusType qcStatusType = QcStatusType.valueOf(row.get(COL_QC_STATUS_TYPE));
                ComparisonOperator comparisonOperator = ComparisonOperator.fromString(row.get(COL_COMPARISON_OPERATOR));
                double threshold = row.getDouble(COL_THRESHOLD);

                ThresholdKey key = new ThresholdKey(sampleType, featureType, featureName, qcStatusType);

                if(!thresholds.containsKey(key))
                {
                    hasUnsupportedThresholds = true;
                    QC_LOGGER.error("  Unsupported threshold: {}", key.toString());
                }

                if(hasUnsupportedThresholds)
                {
                    continue;
                }

                if(seenKeys.contains(key))
                {
                    QC_LOGGER.warn("  Skipping duplicate threshold: {}", key.toString());
                    continue;
                }

                QC_LOGGER.debug("  Overriding threshold: {} default({}) override({})",
                        key.toString(),
                        thresholds.getThreshold(key, false).threshold(),
                        threshold);

                QcThreshold thresholdOverride = QcThreshold.builder(key)
                        .comparisonOperator(comparisonOperator)
                        .threshold(threshold)
                        .build();

                thresholds.overrideThreshold(key, thresholdOverride);

                seenKeys.add(key);
            }

            if(hasUnsupportedThresholds)
            {
                throw new RuntimeException("Failed to read threshold overrides file - unsupported threshold entries found");
            }

            return thresholds.freeze();
        }
    }

    public static void write(String filename, List<QcThreshold> thresholds)
    {
        List<String> columns = new ArrayList<>();
        columns.add(COL_SAMPLE_TYPE);
        columns.add(COL_FEATURE_TYPE);
        columns.add(COL_FEATURE_NAME);
        columns.add(COL_QC_STATUS_TYPE);
        columns.add(COL_COMPARISON_OPERATOR);
        columns.add(COL_THRESHOLD);
        columns.add(COL_OVERRIDABLE);

        DelimFileWriter.write(filename, columns, thresholds, (threshold, row) -> {
            ThresholdKey key = threshold.key();

            row.set(COL_SAMPLE_TYPE, key.sampleType().name());
            row.set(COL_FEATURE_TYPE, key.featureType().name());
            row.set(COL_FEATURE_NAME, key.featureName());
            row.set(COL_QC_STATUS_TYPE, key.qcStatusType().name());

            row.set(COL_COMPARISON_OPERATOR, threshold.determinedElsewhere() ? "NA" : threshold.operator().name());

            row.set(COL_THRESHOLD, threshold.threshold());
            row.set(COL_OVERRIDABLE, threshold.determinedElsewhere() ? Boolean.FALSE.toString() : Boolean.TRUE.toString());
        });
    }
}
