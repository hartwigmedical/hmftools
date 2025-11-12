package com.hartwig.hmftools.qsee.cohort;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.qsee.common.QseeConstants.QC_LOGGER;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.COHORT_FILE_ID;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.COL_FEATURE_NAME;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.COL_FEATURE_TYPE;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.COL_SAMPLE_TYPE;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.COL_SOURCE_TOOL;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.QSEE_FILE_ID;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.utils.file.FileWriterUtils;
import com.hartwig.hmftools.qsee.common.QseeFileCommon;
import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.FeatureKey;
import com.hartwig.hmftools.qsee.feature.FeatureType;
import com.hartwig.hmftools.qsee.feature.SourceTool;

public class CohortPercentilesFile
{
    public static final String COHORT_PERCENTILES_FILE_SUFFIX = COHORT_FILE_ID + "." + QSEE_FILE_ID + ".percentiles.tsv.gz";

    public static final String COL_PERCENTILE_PREFIX = "Pct";
    public static final String COL_PERCENTILE_DELIM = "_";

    public static String generateFilename(final String basePath)
    {
        return checkAddDirSeparator(basePath) + COHORT_PERCENTILES_FILE_SUFFIX;
    }

    public static List<String> getPercentileNames(double[] percentiles)
    {
        return Arrays.stream(percentiles).mapToObj(QseeFileCommon.DECIMAL_FORMAT::format).toList();
    }

    private static List<String> toLines(double[] percentiles, List<FeaturePercentiles> cohortPercentiles)
    {
        List<String> lines = new ArrayList<>();

        StringJoiner header = new StringJoiner(TSV_DELIM);

        header.add(COL_SAMPLE_TYPE);
        header.add(COL_SOURCE_TOOL);
        header.add(COL_FEATURE_TYPE);
        header.add(COL_FEATURE_NAME);

        List<String> percentileNames = getPercentileNames(percentiles);
        for( String percentileName : percentileNames)
        {
            header.add(CohortPercentilesFile.COL_PERCENTILE_PREFIX + "_" + percentileName);
        }

        lines.add(header.toString());

        for(FeaturePercentiles featurePercentiles : cohortPercentiles)
        {
            FeatureKey featureKey = featurePercentiles.featureKey();

            FeatureType featureType = featureKey.type();
            SourceTool sourceTool = featureKey.sourceTool();
            SampleType sampleType = featurePercentiles.sampleType();

            StringJoiner line = new StringJoiner(TSV_DELIM);
            line.add(sampleType.toString());
            line.add(sourceTool.toString());
            line.add(featureType.toString());
            line.add(featureKey.name());

            double[] refValues = featurePercentiles.refValues();

            String refValuesStr = Arrays.stream(refValues)
                    .mapToObj(QseeFileCommon.DECIMAL_FORMAT::format)
                    .collect(Collectors.joining(TSV_DELIM));

            line.add(refValuesStr);

            lines.add(line.toString());
        }

        return lines;
    }

    public static void write(String filename, double[] percentiles, List<FeaturePercentiles> cohortPercentiles) throws IOException
    {
        BufferedWriter writer = FileWriterUtils.createBufferedWriter(filename);

        List<String> lines = toLines(percentiles, cohortPercentiles);
        for(String line : lines)
        {
            writer.write(line);
            writer.newLine();
        }

        writer.close();
    }

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
                    String percentileStr = columnName.split(COL_PERCENTILE_DELIM)[1];
                    percentileList.add(Double.parseDouble(percentileStr));
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
            System.exit(1);
            return null;
        }
    }
}
