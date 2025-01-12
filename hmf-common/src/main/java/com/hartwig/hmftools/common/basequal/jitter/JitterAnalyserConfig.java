package com.hartwig.hmftools.common.basequal.jitter;

import static com.hartwig.hmftools.common.basequal.jitter.JitterAnalyserConstants.DEFAULT_MAX_SINGLE_SITE_ALT_CONTRIBUTION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.region.SpecificRegions.SPECIFIC_REGIONS;
import static com.hartwig.hmftools.common.region.SpecificRegions.loadSpecificRegions;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkCreateOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;

import java.util.List;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.apache.commons.cli.ParseException;

public class JitterAnalyserConfig
{
    public final String SampleId;
    public final RefGenomeVersion RefGenVersion;
    public final String RefGenomeFile;

    public final String RefGenomeMsiFile;
    public final String OutputDir;

    public final int MinMappingQuality;

    public final int MaxSitesPerType;
    public final double MaxSingleSiteAltContribution;

    public final boolean WritePlots;
    public final boolean WriteSiteFile;

    public final List<ChrBaseRegion> SpecificRegions;

    public static final String JITTER_MSI_SITES_FILE = "ref_genome_msi_file";
    public static final String JITTER_MSI_SITES_FILE_DESC = "Path to ref genome MSI sites file";

    public static final String JITTER_MAX_SITES_PER_TYPE = "max_sites_per_type";
    public static final String JITTER_MAX_SITES_PER_TYPE_DESC = "Max number of sites per microsatellite unit / length type";

    private static final String MIN_MAP_QUALITY = "min_map_quality";
    private static final String MAX_SINGLE_SITE_ALT_CONTRIBUTION = "max_site_alt_contribution";

    public static final String JITTER_WRITE_SITE_FILE = "write_msi_site_file";
    public static final String JITTER_WRITE_SITE_FILE_DESC = "Write MSI site file (useful for debugging)";

    public static final int DEFAULT_MIN_MAPPING_QUALITY = 50;
    public static final int DEFAULT_NUM_SITES_PER_TYPE = 5_000;

    public JitterAnalyserConfig(final ConfigBuilder configBuilder) throws ParseException
    {
        SampleId = configBuilder.getValue(SAMPLE);
        RefGenVersion = RefGenomeVersion.from(configBuilder);
        RefGenomeFile = configBuilder.getValue(REF_GENOME);
        RefGenomeMsiFile = configBuilder.getValue(JITTER_MSI_SITES_FILE);
        OutputDir = parseOutputDir(configBuilder);
        MinMappingQuality = configBuilder.getInteger(MIN_MAP_QUALITY);
        MaxSitesPerType = configBuilder.getInteger(JITTER_MAX_SITES_PER_TYPE);
        MaxSingleSiteAltContribution = configBuilder.getDecimal(MAX_SINGLE_SITE_ALT_CONTRIBUTION);
        WritePlots = true;
        WriteSiteFile = configBuilder.hasFlag(JITTER_WRITE_SITE_FILE);
        SpecificRegions = loadSpecificRegions(configBuilder.getValue(SPECIFIC_REGIONS));
    }

    public JitterAnalyserConfig(
            final String sampleId, final RefGenomeVersion refGenVersion, final String refGenomeFile,
            final String refGenomeMsiFile, final String outputDir, double maxSingleSiteAltContribution, boolean writeSiteFile)
    {
        SampleId = sampleId;
        RefGenVersion = refGenVersion;
        RefGenomeFile = refGenomeFile;
        RefGenomeMsiFile = refGenomeMsiFile;
        OutputDir = outputDir;
        MinMappingQuality = JitterAnalyserConfig.DEFAULT_MIN_MAPPING_QUALITY;
        MaxSitesPerType = DEFAULT_NUM_SITES_PER_TYPE;
        MaxSingleSiteAltContribution = maxSingleSiteAltContribution;
        WritePlots = false;
        WriteSiteFile = writeSiteFile;
        SpecificRegions = null;
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(SAMPLE, true, SAMPLE_DESC);

        addRefGenomeConfig(configBuilder, true);

        configBuilder.addPath(JITTER_MSI_SITES_FILE, true, JITTER_MSI_SITES_FILE_DESC);

        addOutputDir(configBuilder);

        configBuilder.addInteger(
                MIN_MAP_QUALITY, "Minimum mapping quality for an alignment to be used", DEFAULT_MIN_MAPPING_QUALITY);

        configBuilder.addInteger(JITTER_MAX_SITES_PER_TYPE, JITTER_MAX_SITES_PER_TYPE_DESC, DEFAULT_NUM_SITES_PER_TYPE);

        configBuilder.addDecimal(
                MAX_SINGLE_SITE_ALT_CONTRIBUTION, "Max percentage a single alt site can contribute",
                DEFAULT_MAX_SINGLE_SITE_ALT_CONTRIBUTION);

        configBuilder.addFlag(JITTER_WRITE_SITE_FILE, JITTER_WRITE_SITE_FILE_DESC);
    }

    public boolean isValid()
    {
        checkCreateOutputDir(OutputDir);
        return true;
    }
}
