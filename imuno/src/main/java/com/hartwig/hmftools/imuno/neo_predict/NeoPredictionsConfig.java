package com.hartwig.hmftools.imuno.neo_predict;

import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.DELIMITER;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.imuno.common.ImunoCommon.IM_LOGGER;
import static com.hartwig.hmftools.imuno.common.ImunoCommon.LOG_DEBUG;
import static com.hartwig.hmftools.imuno.common.ImunoCommon.loadSampleDataFile;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.imuno.neo.SampleData;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class NeoPredictionsConfig
{
    public final List<SampleData> Samples;

    public final boolean WriteCohortFile;
    public final String SampleDataDir;
    public final String OutputDir;

    public static final String SAMPLE = "sample";
    public static final String SAMPLE_DATA_DIR = "sample_data_dir";
    public static final String WRITE_COHORT_FILE = "write_cohort_file";

    public static final String PREDICTIONS_FILE_ID = ".mcf.predictions.csv";

    public NeoPredictionsConfig(final CommandLine cmd)
    {
        Samples = Lists.newArrayList();

        final String sampleIdConfig = cmd.getOptionValue(SAMPLE);

        if(sampleIdConfig.contains(".csv"))
        {
            loadSampleDataFile(sampleIdConfig, Samples);
        }
        else
        {
            if(sampleIdConfig.contains(DELIMITER))
            {
                SampleData sample = SampleData.fromCsv(sampleIdConfig);

                if(sample == null)
                {
                    IM_LOGGER.error("invalid sample data: {}", sampleIdConfig);
                }
                else
                {
                    Samples.add(sample);
                }
            }
            else
            {
                Samples.add(new SampleData(sampleIdConfig, "", Lists.newArrayList()));
            }
        }

        SampleDataDir = cmd.getOptionValue(SAMPLE_DATA_DIR);
        OutputDir = parseOutputDir(cmd);

        WriteCohortFile = cmd.hasOption(WRITE_COHORT_FILE);
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(SAMPLE, true, "Sample - Id(s) separated by ';' or CSV file");
        options.addOption(SAMPLE_DATA_DIR, true, "SV fusion file (single sample or cohort)");
        options.addOption(WRITE_COHORT_FILE, false, "Write cohort files for multiple samples");
        options.addOption(OUTPUT_DIR, true, "Output directory");
        options.addOption(LOG_DEBUG, false, "Log verbose");
        DatabaseAccess.addDatabaseCmdLineArgs(options);
    }
}
