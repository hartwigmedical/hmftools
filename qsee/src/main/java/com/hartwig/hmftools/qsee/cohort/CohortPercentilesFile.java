package com.hartwig.hmftools.qsee.cohort;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.COL_FEATURE_NAME;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.COL_FEATURE_TYPE;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.COL_SAMPLE_TYPE;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.COL_SOURCE_TOOL;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.QSEE_FILE_ID;
import static com.hartwig.hmftools.qsee.common.QseeConstants.QC_LOGGER;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.utils.file.FileWriterUtils;
import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.FeatureKey;
import com.hartwig.hmftools.qsee.feature.FeatureType;
import com.hartwig.hmftools.qsee.feature.SourceTool;

public class CohortPercentilesFile
{
    public static final String COHORT_PERCENTILES_FILE_SUFFIX = "cohort." + QSEE_FILE_ID + ".percentiles.tsv.gz";

    public static final String COL_PERCENTILE_PREFIX = "Pct_";

    public static final DecimalFormat PERCENTILE_FORMAT = new DecimalFormat("0.########");
    public static final DecimalFormat REF_VALUE_FORMAT  = new DecimalFormat("0.########");

    public static CohortPercentiles read(String filename)
    {
        try
        {
            List<String> lines = FileWriterUtils.readLines(filename);

            String header = lines.get(0);
            Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);
            lines.remove(0);

            int indexSampleType = fieldsIndexMap.get(COL_SAMPLE_TYPE);
            int indexFeatureType = fieldsIndexMap.get(COL_FEATURE_TYPE);
            int indexFeatureName = fieldsIndexMap.get(COL_FEATURE_NAME);
            int indexSourceTool = fieldsIndexMap.get(COL_SOURCE_TOOL);

            List<Integer> indexesPercentile = new ArrayList<>();
            List<Double> percentileList = new ArrayList<>();
            for(Map.Entry<String, Integer> entry : fieldsIndexMap.entrySet())
            {
                String columnName = entry.getKey();
                if(columnName.startsWith(COL_PERCENTILE_PREFIX))
                {
                    indexesPercentile.add(entry.getValue());
                    percentileList.add(Double.parseDouble(columnName.substring(COL_PERCENTILE_PREFIX.length())));
                }
            }

            double[] percentiles = percentileList.stream().mapToDouble(Double::doubleValue).toArray();

            CohortPercentiles cohortPercentiles = new CohortPercentiles();

            for(String line : lines)
            {
                String[] fields = line.split(TSV_DELIM);
                SampleType sampleType = SampleType.valueOf(fields[indexSampleType]);

                String featureName = fields[indexFeatureName];
                FeatureType featureType = FeatureType.valueOf(fields[indexFeatureType]);
                SourceTool sourceTool = SourceTool.valueOf(fields[indexSourceTool]);
                FeatureKey featureKey = new FeatureKey(featureName, featureType, sourceTool);

                double[] refValues = indexesPercentile.stream()
                        .mapToDouble(i -> Double.parseDouble(fields[i]))
                        .toArray();

                FeaturePercentiles featurePercentiles = new FeaturePercentiles(sampleType, featureKey, percentiles, refValues);

                cohortPercentiles.add(featurePercentiles);
            }

            return cohortPercentiles;
        }
        catch(IOException e)
        {
            QC_LOGGER.error("Failed to load cohort percentiles file: {}", filename, e);
            return null;
        }
    }
}
