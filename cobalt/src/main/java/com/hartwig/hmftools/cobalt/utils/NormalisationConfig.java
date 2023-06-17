package com.hartwig.hmftools.cobalt.utils;

import static com.hartwig.hmftools.common.genome.gc.GCProfileFactory.GC_PROFILE;
import static com.hartwig.hmftools.common.genome.gc.GCProfileFactory.addGcProfilePath;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION_CFG_DESC;
import static com.hartwig.hmftools.common.utils.ConfigUtils.SAMPLE_ID_FILE;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addSampleIdFile;
import static com.hartwig.hmftools.common.utils.ConfigUtils.loadSampleIdsFile;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;

import java.util.List;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class NormalisationConfig
{
    public final List<String> SampleIds;
    public final String AmberDir;
    public final String CobaltDir;
    public final String TargetRegionsBed;
    public final String GcProfile;
    public final String OutputDir;
    public final String OutputId;
    public final RefGenomeVersion RefGenVersion;

    private static final String COBALT_DIR = "cobalt_dir";
    private static final String AMBER_DIR = "amber_dir";
    private static final String TARGET_REGIONS_BED = "target_regions_bed";

    public NormalisationConfig(final CommandLine cmd)
    {
        SampleIds = loadSampleIdsFile(cmd);
        CobaltDir = cmd.getOptionValue(COBALT_DIR);
        AmberDir = cmd.getOptionValue(AMBER_DIR);
        GcProfile = cmd.getOptionValue(GC_PROFILE);
        TargetRegionsBed = cmd.getOptionValue(TARGET_REGIONS_BED);
        RefGenVersion = RefGenomeVersion.from(cmd);
        OutputDir = parseOutputDir(cmd);
        OutputId = cmd.getOptionValue(OUTPUT_ID);
    }

    public static void addCommandLineOptions(final Options options)
    {
        options.addOption(SAMPLE_ID_FILE, true, "Sample IDs file");
        options.addOption(AMBER_DIR, true, "Path to amber files");
        options.addOption(COBALT_DIR, true, "Path to cobalt files");
        options.addOption(REF_GENOME_VERSION, true, REF_GENOME_VERSION_CFG_DESC);
        options.addOption(TARGET_REGIONS_BED, true, "Target regions BED file");
        addGcProfilePath(options);
        addSampleIdFile(options);
        addOutputOptions(options);
        addLoggingOptions(options);
    }
}
