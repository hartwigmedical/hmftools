package com.hartwig.hmftools.svprep.append;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.region.SpecificRegions.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.region.SpecificRegions.loadSpecificRegions;
import static com.hartwig.hmftools.common.samtools.BamUtils.addValidationStringencyOption;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LOG_READ_IDS;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LOG_READ_IDS_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.parseLogReadIds;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;
import static com.hartwig.hmftools.svprep.SvConfig.BAM_FILE;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.svprep.reads.ReadFilterConfig;
import com.hartwig.hmftools.svprep.reads.ReadFilters;

import org.apache.commons.cli.ParseException;

public class AppendConfig
{
    public final String SampleId;
    public final String BamFile;
    public final String InputVcf;
    public final String OutputVcf;

    public final String RefGenomeFile;
    public final RefGenomeVersion RefGenVersion;

    public final ReadFilters ReadFiltering;

    public final int Threads;

    // debug
    public final List<String> LogReadIds;
    public final List<ChrBaseRegion> SpecificRegions;

    public static final String INPUT_VCF = "input_vcf";
    public static final String OUTPUT_VCF = "output_vcf";

    public AppendConfig(final ConfigBuilder configBuilder)
    {
        SampleId = configBuilder.getValue(SAMPLE);
        BamFile = configBuilder.getValue(BAM_FILE);
        RefGenomeFile = configBuilder.getValue(REF_GENOME);
        InputVcf = configBuilder.getValue(INPUT_VCF);
        OutputVcf = configBuilder.getValue(OUTPUT_VCF);

        RefGenVersion = RefGenomeVersion.from(configBuilder);
        SV_LOGGER.info("refGenome({}), bam({})", RefGenVersion, BamFile);

        ReadFiltering = new ReadFilters(ReadFilterConfig.from(configBuilder));

        SpecificRegions = Lists.newArrayList();

        try
        {
            SpecificRegions.addAll(loadSpecificRegions(configBuilder));
        }
        catch(ParseException e)
        {
            System.exit(1);
        }

        LogReadIds = parseLogReadIds(configBuilder);

        Threads = parseThreads(configBuilder);
    }

    public boolean isValid()
    {
        return checkFileExists(BamFile) && checkFileExists(InputVcf) && checkFileExists(RefGenomeFile);
    }

    private static boolean checkFileExists(final String filename)
    {
        if(!Files.exists(Paths.get(filename)))
        {
            SV_LOGGER.error("missing required file({})", filename);
            return false;
        }

        return true;
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(SAMPLE, true, SAMPLE_DESC);
        configBuilder.addPath(BAM_FILE, true, "BAM file location");
        configBuilder.addPath(INPUT_VCF, true, "Input VCF");
        configBuilder.addConfigItem(OUTPUT_VCF, true, "Output VCF");
        addRefGenomeConfig(configBuilder, true);
        addValidationStringencyOption(configBuilder);
        ReadFilterConfig.addConfig(configBuilder);

        addSpecificChromosomesRegionsConfig(configBuilder);
        configBuilder.addConfigItem(LOG_READ_IDS, false, LOG_READ_IDS_DESC);
        addOutputOptions(configBuilder);
        ConfigUtils.addLoggingOptions(configBuilder);
        addThreadOptions(configBuilder);
    }
}
