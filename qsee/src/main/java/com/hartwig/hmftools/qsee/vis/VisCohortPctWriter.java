package com.hartwig.hmftools.qsee.vis;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.qsee.common.QseeConstants.APP_NAME;
import static com.hartwig.hmftools.qsee.common.QseeConstants.QC_LOGGER;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.COL_FEATURE_NAME;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.COL_FEATURE_TYPE;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.COL_SAMPLE_TYPE;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.COL_SOURCE_TOOL;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.QSEE_FILE_ID;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.qsee.cohort.CohortPercentiles;
import com.hartwig.hmftools.qsee.cohort.CohortPercentilesFile;
import com.hartwig.hmftools.qsee.cohort.FeaturePercentiles;
import com.hartwig.hmftools.qsee.cohort.NamedPercentile;
import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.FeatureKey;

public class VisCohortPctWriter
{
    private final VisPrepConfig mConfig;

    public VisCohortPctWriter(VisPrepConfig config)
    {
        mConfig = config;
    }

    private static List<Map<NamedPercentile, Double>> getNamedPercentiles(CohortPercentiles cohortPercentiles)
    {
        List<Map<NamedPercentile, Double>> namedPercentilesList = new ArrayList<>();
        for(SampleType sampleType : cohortPercentiles.getData().keySet())
        {
            Map<FeatureKey, FeaturePercentiles> featurePercentilesMap = cohortPercentiles.getData().get(sampleType);
            for(FeatureKey featureKey : featurePercentilesMap.keySet())
            {
                FeaturePercentiles featurePercentiles = featurePercentilesMap.get(featureKey);
                Map<NamedPercentile, Double> namedPercentiles = featurePercentiles.getNamedPercentileValues();
                namedPercentilesList.add(namedPercentiles);
            }
        }

        return namedPercentilesList;
    }

    private void writeToFile(String outputFile, List<Map<NamedPercentile, Double>> namedPercentilesList)
    {
        try(BufferedWriter writer = createBufferedWriter(outputFile))
        {
            QC_LOGGER.info("Writing cohort named percentiles data to: {}", outputFile);

            StringJoiner header = new StringJoiner(TSV_DELIM);
            header.add(COL_SAMPLE_TYPE);
            header.add(COL_SOURCE_TOOL);
            header.add(COL_FEATURE_TYPE);
            header.add(COL_FEATURE_NAME);

            for(NamedPercentile namedPercentile : NamedPercentile.values())
            {
                header.add(namedPercentile.name());
            }

            writer.write(header.toString());
            writer.newLine();

            for(Map<NamedPercentile, Double> pctFeatureValueMap : namedPercentilesList)
            {
                StringJoiner line = new StringJoiner(TSV_DELIM);
                for(NamedPercentile namedPercentile : pctFeatureValueMap.keySet())
                {
                    String featureValue = CohortPercentilesFile.PERCENTILE_FORMAT.format(pctFeatureValueMap.get(namedPercentile));
                    line.add(featureValue);
                }
                writer.write(line.toString());
                writer.newLine();
            }
        }
        catch(IOException e)
        {
            QC_LOGGER.error("Failed to write to file: {}", outputFile, e);
            System.exit(1);
        }
    }

    public void run()
    {
        CohortPercentiles cohortPercentiles = CohortPercentilesFile.read(mConfig.CohortPercentilesFile);

        List<Map<NamedPercentile, Double>> namedPercentilesList = getNamedPercentiles(cohortPercentiles);

        String outputFile = checkAddDirSeparator(mConfig.CommonPrep.OutputDir) + "cohort." + QSEE_FILE_ID + ".vis.named_percentiles.tsv.gz";
        writeToFile(outputFile, namedPercentilesList);
    }

    public static void main(String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        VisPrepConfig.registerConfig(configBuilder);
        configBuilder.checkAndParseCommandLine(args);

        VisPrepConfig config = new VisPrepConfig(configBuilder);
        new VisCohortPctWriter(config).run();
    }
}
