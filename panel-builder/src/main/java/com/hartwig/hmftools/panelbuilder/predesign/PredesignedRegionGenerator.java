package com.hartwig.hmftools.panelbuilder.predesign;

import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR_DESC;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkCreateOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

// Generates Hartwig predesigned regions which can be used as "custom region" inputs in PanelBuilder.
public class PredesignedRegionGenerator
{
    private final String mMsiSitesFile;
    private final String mOutputDir;

    private static final String CFG_MSI_SITES_FILE = "msi_sites";
    private static final String DESC_MSI_SITES_FILE = "Microsatellite instability positions TSV file";

    private PredesignedRegionGenerator(final ConfigBuilder configBuilder)
    {
        mMsiSitesFile = configBuilder.getValue(CFG_MSI_SITES_FILE);
        mOutputDir = parseOutputDir(configBuilder);
    }

    private void run()
    {
        checkCreateOutputDir(mOutputDir);

        MsiSites.generateRegions(mMsiSitesFile, mOutputDir);
    }

    public static void main(final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder();

        configBuilder.addPath(CFG_MSI_SITES_FILE, true, DESC_MSI_SITES_FILE);
        configBuilder.addConfigItem(OUTPUT_DIR, true, OUTPUT_DIR_DESC);
        addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        PredesignedRegionGenerator predesignedRegionGenerator = new PredesignedRegionGenerator(configBuilder);
        predesignedRegionGenerator.run();
    }
}
