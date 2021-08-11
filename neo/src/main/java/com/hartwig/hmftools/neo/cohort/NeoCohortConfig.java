package com.hartwig.hmftools.neo.cohort;

import static com.hartwig.hmftools.common.utils.ConfigUtils.LOG_DEBUG;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.neo.NeoCommon.loadSampleIdsFile;
import static com.hartwig.hmftools.neo.bind.BinderConfig.OUTPUT_ID;
import static com.hartwig.hmftools.neo.cohort.CohortWriteDetail.NEOEPITOPE;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.ConfigUtils;
import com.hartwig.hmftools.neo.bind.BinderConfig;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class NeoCohortConfig
{
    public final String OutputDir;
    public final String OutputId;
    public final String NeoDataDir;
    public final String LilacDataDir;
    public final String McfPredictionsDir;
    public final List<String> SampleIds;

    public final CohortWriteDetail WriteDetail;

    public final double McfSumFactor;
    public final int Threads;

    public static final String SAMPLE_ID_FILE = "sample_id_file";
    public static final String NEO_DATA_DIR = "neo_data_dir";
    public static final String LILAC_DATA_DIR = "lilac_data_dir";
    public static final String PREDICTION_DATA_DIR = "mcf_prediction_dir";
    public static final String THREADS = "threads";

    private static final String MCF_SUM_FACTOR = "mcf_sum_factor";

    private static final String WRITE_DETAIL = "write_detail";

    public NeoCohortConfig(final CommandLine cmd)
    {
        SampleIds = Lists.newArrayList();
        loadSampleIdsFile(cmd.getOptionValue(SAMPLE_ID_FILE), SampleIds);

        NeoDataDir = cmd.getOptionValue(NEO_DATA_DIR);
        McfPredictionsDir = cmd.getOptionValue(PREDICTION_DATA_DIR);
        LilacDataDir = cmd.getOptionValue(LILAC_DATA_DIR);

        McfSumFactor = Double.parseDouble(cmd.getOptionValue(MCF_SUM_FACTOR, "2"));

        OutputDir = parseOutputDir(cmd);
        OutputId = cmd.getOptionValue(OUTPUT_ID);
        Threads = Integer.parseInt(cmd.getOptionValue(THREADS, "0"));
        WriteDetail = CohortWriteDetail.valueOf(cmd.getOptionValue(WRITE_DETAIL, String.valueOf(NEOEPITOPE)));
    }

    public String formFilename(final String fileId)
    {
        if(OutputId == null || OutputId.isEmpty())
            return String.format("%sneo_cohort_%s.csv", OutputDir, fileId);
        else
            return String.format("%sneo_cohort_%s_%s.csv", OutputDir, OutputId, fileId);
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(SAMPLE_ID_FILE, true, "SampleId file");
        options.addOption(NEO_DATA_DIR, true, "Directory for sample neo-epitope files");
        options.addOption(PREDICTION_DATA_DIR, true, "Directory for sample prediction result files");
        options.addOption(LILAC_DATA_DIR, true, "Directory for Lilac coverage files");

        BinderConfig.addCmdLineArgs(options);
        ConfigUtils.addLoggingOptions(options);
        options.addOption(MCF_SUM_FACTOR, true, "Affinity sum factor");
        options.addOption(THREADS, true, "Thread count");

        options.addOption(WRITE_DETAIL, false, "Write all peptide scores");
        options.addOption(OUTPUT_DIR, true, "Output directory");
        options.addOption(OUTPUT_ID, true, "Output file ID");
    }

}
