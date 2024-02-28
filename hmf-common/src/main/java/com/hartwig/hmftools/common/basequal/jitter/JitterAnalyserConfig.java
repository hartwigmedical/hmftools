package com.hartwig.hmftools.common.basequal.jitter;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeVersion;
import static com.hartwig.hmftools.common.region.SpecificRegions.SPECIFIC_REGIONS;
import static com.hartwig.hmftools.common.region.SpecificRegions.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.region.SpecificRegions.loadSpecificRegions;
import static com.hartwig.hmftools.common.samtools.BamUtils.addValidationStringencyOption;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkCreateOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;

import java.util.List;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.samtools.BamUtils;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.apache.commons.cli.ParseException;

import htsjdk.samtools.ValidationStringency;

public class JitterAnalyserConfig
{
    public final String SampleId;
    public final String BamPath;
    public final String RefGenomeFile;

    public final String RefGenomeMicrosatelliteFile;
    public final String OutputDir;

    public final int MinMappingQuality;

    public final int MaxSitesPerType;

    public final ValidationStringency BamStringency;
    public final int Threads;

    public final List<ChrBaseRegion> SpecificRegions;

    private static final String REF_GENOME_MICROSATELLITES = "ref_genome_microsatellites";
    private static final String MIN_MAP_QUALITY = "min_map_quality";
    private static final String MAX_SITES_PER_TYPE = "max_sites_per_type";
    public static final int DEFAULT_MIN_MAPPING_QUALITY = 50;
    public static final int DEFAULT_NUM_SITES_PER_TYPE = 5000;

    public JitterAnalyserConfig(final ConfigBuilder configBuilder) throws ParseException
    {
        SampleId = configBuilder.getValue(SAMPLE);
        BamPath = configBuilder.getValue("bam");
        RefGenomeFile = configBuilder.getValue(REF_GENOME);
        RefGenomeMicrosatelliteFile = configBuilder.getValue(REF_GENOME_MICROSATELLITES);
        OutputDir = parseOutputDir(configBuilder);
        Threads = parseThreads(configBuilder);
        MinMappingQuality = configBuilder.getInteger(MIN_MAP_QUALITY);
        MaxSitesPerType = configBuilder.getInteger(MAX_SITES_PER_TYPE);
        BamStringency = BamUtils.validationStringency(configBuilder);
        SpecificRegions = loadSpecificRegions(configBuilder.getValue(SPECIFIC_REGIONS));
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(SAMPLE, true, "sample id");
        configBuilder.addPath("bam", true, "path to bam file");

        addRefGenomeVersion(configBuilder);
        configBuilder.addPath(REF_GENOME, true, REF_GENOME_CFG_DESC + ", required when using CRAM files");

        configBuilder.addPath(REF_GENOME_MICROSATELLITES, true, "path to ref genome microsatellites tsv");

        addOutputDir(configBuilder);

        configBuilder.addInteger(
                MIN_MAP_QUALITY, "Minimum mapping quality for an alignment to be used", DEFAULT_MIN_MAPPING_QUALITY);

        configBuilder.addInteger(MAX_SITES_PER_TYPE, "max number of sites per microsatellite unit / length type", DEFAULT_NUM_SITES_PER_TYPE);

        addThreadOptions(configBuilder);
        addValidationStringencyOption(configBuilder);
        addLoggingOptions(configBuilder);

        addSpecificChromosomesRegionsConfig(configBuilder);
    }

    public boolean isValid()
    {
        checkCreateOutputDir(OutputDir);
        return true;
    }
}
