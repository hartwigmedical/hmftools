package com.hartwig.hmftools.qsee.cohort;

import static java.util.Objects.requireNonNull;

import static com.hartwig.hmftools.qsee.common.QseeConstants.QC_LOGGER;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.COL_FEATURE_NAME;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.COL_FEATURE_TYPE;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.COL_SAMPLE_TYPE;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.COL_SOURCE_TOOL;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.hartwig.hmftools.common.utils.file.DelimFileReader;
import com.hartwig.hmftools.common.utils.file.DelimFileWriter;
import com.hartwig.hmftools.qsee.common.QseeFileCommon;
import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.FeatureKey;
import com.hartwig.hmftools.qsee.feature.FeatureType;
import com.hartwig.hmftools.qsee.feature.SourceTool;

import org.jetbrains.annotations.Nullable;

public class CohortPercentilesFile
{
    public static final String COHORT_PERCENTILES_FILE_CFG = "cohort_percentiles_file";
    public static final String COHORT_PERCENTILES_FILE_CFG_DESC = "Path to the cohort percentiles file";

    public static final String COL_PERCENTILE_PREFIX = "Pct";
    public static final String COL_PERCENTILE_DELIM = "_";

    public static String generateFilename(String basePath, @Nullable String outputId)
    {
        return QseeFileCommon.generateCohortFilename(basePath, "percentiles", outputId, "tsv.gz");
    }

    public static String generateFilename(String basePath)
    {
        return generateFilename(basePath, null);
    }

    public static List<String> getPercentileNames(double[] percentiles)
    {
        return Arrays.stream(percentiles).mapToObj(QseeFileCommon.DECIMAL_FORMAT::format).toList();
    }

    public static void write(String filename, double[] percentiles, List<FeaturePercentiles> cohortPercentiles) throws IOException
    {
        List<String> percentileNames = getPercentileNames(percentiles);

        List<String> columns = new ArrayList<>();
        columns.add(COL_SAMPLE_TYPE);
        columns.add(COL_SOURCE_TOOL);
        columns.add(COL_FEATURE_TYPE);
        columns.add(COL_FEATURE_NAME);

        for(String percentileName : percentileNames)
        {
            columns.add(COL_PERCENTILE_PREFIX + COL_PERCENTILE_DELIM + percentileName);
        }

        DelimFileWriter.write(filename, columns, cohortPercentiles, (featurePercentiles, row) ->
        {
            FeatureKey featureKey = featurePercentiles.featureKey();

            row.set(COL_SAMPLE_TYPE, featurePercentiles.sampleType().name());
            row.set(COL_SOURCE_TOOL, featureKey.sourceTool().name());
            row.set(COL_FEATURE_TYPE, featureKey.type().name());
            row.set(COL_FEATURE_NAME, featureKey.name());

            double[] refValues = featurePercentiles.refValues();
            for(int i = 0; i < percentileNames.size(); i++)
            {
                String col = COL_PERCENTILE_PREFIX + COL_PERCENTILE_DELIM + percentileNames.get(i);
                row.set(col, refValues[i], QseeFileCommon.DECIMAL_FORMAT);
            }
        });
    }

    public static CohortPercentiles read(String filename)
    {
        try(DelimFileReader reader = new DelimFileReader(filename))
        {
            int indexSampleType = requireNonNull(reader.getColumnIndex(COL_SAMPLE_TYPE));
            int indexFeatureType = requireNonNull(reader.getColumnIndex(COL_FEATURE_TYPE));
            int indexFeatureName = requireNonNull(reader.getColumnIndex(COL_FEATURE_NAME));
            int indexSourceTool = requireNonNull(reader.getColumnIndex(COL_SOURCE_TOOL));

            List<String> percentileColumns = new ArrayList<>();
            List<Double> percentileList = new ArrayList<>();

            for(String columnName : reader.getColumnNames())
            {
                if(columnName.startsWith(COL_PERCENTILE_PREFIX))
                {
                    percentileColumns.add(columnName);
                    String percentileStr = columnName.split(COL_PERCENTILE_DELIM, 2)[1];
                    percentileList.add(Double.parseDouble(percentileStr));
                }
            }

            double[] percentiles = percentileList.stream().mapToDouble(Double::doubleValue).toArray();

            CohortPercentiles cohortPercentiles = new CohortPercentiles();

            for(DelimFileReader.Row row : reader)
            {
                SampleType sampleType = SampleType.valueOf(row.get(indexSampleType));

                String featureName = row.get(indexFeatureName);
                FeatureType featureType = FeatureType.valueOf(row.get(indexFeatureType));
                SourceTool sourceTool = SourceTool.valueOf(row.get(indexSourceTool));
                FeatureKey featureKey = new FeatureKey(featureName, featureType, sourceTool);

                double[] percentileValues = new double[percentileColumns.size()];
                for(int i = 0; i < percentileColumns.size(); i++)
                {
                    percentileValues[i] = row.getDouble(percentileColumns.get(i));
                }

                FeaturePercentiles featurePercentiles = new FeaturePercentiles(sampleType, featureKey, percentiles, percentileValues);
                cohortPercentiles.add(featurePercentiles);
            }

            return cohortPercentiles;
        }
        catch(Exception e)
        {
            QC_LOGGER.error("Failed to load cohort percentiles file: {}", filename, e);
            System.exit(1);
            return null;
        }
    }
}
