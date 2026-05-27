package com.hartwig.hmftools.amber.utils;

import static com.hartwig.hmftools.amber.AmberConfig.AMB_LOGGER;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;

import java.util.List;

import com.hartwig.hmftools.amber.contamination.TumorContamination;
import com.hartwig.hmftools.amber.contamination.TumorContaminationFile;
import com.hartwig.hmftools.amber.contamination.TumorContaminationModel;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.jetbrains.annotations.NotNull;

public class ContaminationModelTester
{
    private static final String CONTAMINATION_FILE = "contamination_file";
    private static final String HET_SITES_COUNT = "het_sites_count";

    private final String mContaminationFile;
    private final String mOutputDir;
    private final int mHetSitesCount;

    public ContaminationModelTester(final ConfigBuilder configBuilder)
    {
        mContaminationFile = configBuilder.getValue(CONTAMINATION_FILE);
        mOutputDir = parseOutputDir(configBuilder);
        mHetSitesCount = configBuilder.getInteger(HET_SITES_COUNT);
    }

    public void run()
    {
        AMB_LOGGER.info("running Amber contamination model");

        List<TumorContamination> tumorContaminations = TumorContaminationFile.read(mContaminationFile);

        TumorContaminationModel model = new TumorContaminationModel();
        model.calcContamination(tumorContaminations, mHetSitesCount);

        AMB_LOGGER.info("Amber model runs complete");
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder();

        configBuilder.addPath(CONTAMINATION_FILE, true, "Input germline locations file");
        configBuilder.addInteger(HET_SITES_COUNT,"Count of het sites", 0);
        addOutputDir(configBuilder);
        addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        ContaminationModelTester modelTester = new ContaminationModelTester(configBuilder);
        modelTester.run();
    }
}
