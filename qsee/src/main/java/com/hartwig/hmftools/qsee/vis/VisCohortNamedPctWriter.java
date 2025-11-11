package com.hartwig.hmftools.qsee.vis;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
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
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.hartwig.hmftools.qsee.cohort.CohortPercentiles;
import com.hartwig.hmftools.qsee.cohort.CohortPercentilesFile;
import com.hartwig.hmftools.qsee.cohort.FeaturePercentiles;
import com.hartwig.hmftools.qsee.cohort.NamedPercentile;
import com.hartwig.hmftools.qsee.common.QseeFileCommon;
import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.FeatureKey;

public class VisCohortNamedPctWriter
{
    private final CohortPercentiles mCohortPercentiles;
    private final String mOutputDir;

    public VisCohortNamedPctWriter(CohortPercentiles cohortPercentiles, String outputDir)
    {
        mCohortPercentiles = cohortPercentiles;
        mOutputDir = outputDir;
    }

    public static String generateFilename(String basePath)
    {
        return checkAddDirSeparator(basePath) + COHORT_FILE_ID + "." + QSEE_FILE_ID + ".vis.named_percentiles.tsv.gz";
    }

    private List<VisFeatureNamedPct> getNamedPercentiles()
    {
        List<VisFeatureNamedPct> cohortNamedPercentiles = new ArrayList<>();
        for(SampleType sampleType : mCohortPercentiles.getData().keySet())
        {
            Map<FeatureKey, FeaturePercentiles> featurePercentilesMap = mCohortPercentiles.getData().get(sampleType);
            for(FeatureKey featureKey : featurePercentilesMap.keySet())
            {
                FeaturePercentiles featurePercentiles = featurePercentilesMap.get(featureKey);
                VisFeatureNamedPct featureNamedPercentiles = new VisFeatureNamedPct(sampleType, featureKey, featurePercentiles);
                cohortNamedPercentiles.add(featureNamedPercentiles);
            }
        }

        return cohortNamedPercentiles;
    }

    private void writeToFile(String outputFile, List<VisFeatureNamedPct> cohortNamedPercentiles)
    {
        try(BufferedWriter writer = createBufferedWriter(outputFile))
        {
            QC_LOGGER.info("Writing cohort named percentiles to: {}", outputFile);

            StringJoiner header = new StringJoiner(TSV_DELIM);
            header.add(COL_SAMPLE_TYPE);
            header.add(COL_SOURCE_TOOL);
            header.add(COL_FEATURE_TYPE);
            header.add(COL_FEATURE_NAME);

            for(NamedPercentile namedPercentile : NamedPercentile.values())
            {
                header.add(CohortPercentilesFile.COL_PERCENTILE_PREFIX + namedPercentile.simpleName());
            }

            writer.write(header.toString());
            writer.newLine();

            for(VisFeatureNamedPct namedPercentiles : cohortNamedPercentiles)
            {
                StringJoiner line = new StringJoiner(TSV_DELIM);

                FeatureKey featureKey = namedPercentiles.featureKey();
                line.add(namedPercentiles.sampleType().toString());
                line.add(featureKey.sourceTool().toString());
                line.add(featureKey.type().toString());
                line.add(featureKey.name());

                for(NamedPercentile namedPercentile : namedPercentiles.pctRefValues().keySet())
                {
                    double featureValue = namedPercentiles.pctRefValues().get(namedPercentile);
                    String featureValueStr = QseeFileCommon.DECIMAL_FORMAT.format(featureValue);
                    line.add(featureValueStr);
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
        List<VisFeatureNamedPct> namedPercentilesList = getNamedPercentiles();
        String outputFile = generateFilename(mOutputDir);
        writeToFile(outputFile, namedPercentilesList);
    }
}
