package com.hartwig.hmftools.redux.jitter;

import static com.hartwig.hmftools.redux.ReduxConfig.SEQUENCING_TYPE;
import static com.hartwig.hmftools.redux.ReduxConstants.BQR_MIN_MAP_QUAL;
import static com.hartwig.hmftools.redux.jitter.JitterConstants.DEFAULT_MAX_SINGLE_SITE_ALT_CONTRIBUTION;
import static com.hartwig.hmftools.common.sequencing.SequencingType.BIOMODAL;
import static com.hartwig.hmftools.common.sequencing.SequencingType.ILLUMINA;
import static com.hartwig.hmftools.common.sequencing.SequencingType.SBX;
import static com.hartwig.hmftools.common.sequencing.SequencingType.ULTIMA;

import java.util.Collections;
import java.util.EnumSet;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.bam.ConsensusType;
import com.hartwig.hmftools.common.region.SpecificRegions;
import com.hartwig.hmftools.common.sequencing.SequencingType;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.jetbrains.annotations.Nullable;

public class MsJitterConfig
{
    public final String SampleId;
    public final RefGenomeVersion RefGenVersion;
    public final String RefGenomeFile;
    public final boolean UsesDuplexUMIs;

    public final String RefGenomeMsiFile;
    public final String OutputDir;

    public final int MinMappingQuality;
    public final int MaxSitesPerType;
    public final double MaxSingleSiteAltContribution;
    public final SpecificRegions SpecificChrRegions;

    public static boolean JITTER_APPLY_BQ_FILTER = true;

    public final boolean WritePlots;
    public final boolean WriteSiteFile;

    public static final String JITTER_WRITE_MSI_PLOTS = "write_msi_plots";
    public static final String JITTER_WRITE_MSI_PLOTS_DESC = "Generate MSI jitter plots";

    public static final String JITTER_MSI_SITES_FILE = "ref_genome_msi_file";
    public static final String JITTER_MSI_SITES_FILE_DESC = "Path to ref genome MSI sites file";

    public static final String JITTER_MAX_SITES_PER_TYPE = "msi_max_sites_per_type";
    public static final String JITTER_MAX_SITES_PER_TYPE_DESC = "Max number of sites per microsatellite unit / length type";

    private static final String NO_BASE_QUAL_FILTER = "msi_no_bq_filter";

    private static final String MIN_MAP_QUALITY = "msi_min_map_quality";
    private static final String MAX_SINGLE_SITE_ALT_CONTRIBUTION = "max_site_alt_contribution";

    public static final String JITTER_WRITE_SITE_FILE = "write_msi_site_file";
    public static final String JITTER_WRITE_SITE_FILE_DESC = "Write MSI site file (useful for debugging)";

    public static final int DEFAULT_MIN_MAPPING_QUALITY = BQR_MIN_MAP_QUAL;
    public static final int DEFAULT_NUM_SITES_PER_TYPE = 5_000;

    private MsJitterConfig(
            final String sampleId, final String refGenomeFile, final RefGenomeVersion refGenVersion,
            boolean usesDuplexUMIs, final String outputDir, final ConfigBuilder configBuilder)
    {
        SampleId = sampleId;
        RefGenomeFile = refGenomeFile;
        RefGenVersion = refGenVersion;
        UsesDuplexUMIs = usesDuplexUMIs;
        OutputDir = outputDir;

        RefGenomeMsiFile = configBuilder.getValue(JITTER_MSI_SITES_FILE);
        MinMappingQuality = configBuilder.getInteger(MIN_MAP_QUALITY);
        MaxSitesPerType = configBuilder.getInteger(JITTER_MAX_SITES_PER_TYPE);
        MaxSingleSiteAltContribution = configBuilder.getDecimal(MAX_SINGLE_SITE_ALT_CONTRIBUTION);

        SpecificChrRegions = SpecificRegions.from(configBuilder, false);

        WritePlots = configBuilder.hasFlag(JITTER_WRITE_MSI_PLOTS);
        WriteSiteFile = configBuilder.hasFlag(JITTER_WRITE_SITE_FILE);

        JITTER_APPLY_BQ_FILTER = !configBuilder.hasFlag(NO_BASE_QUAL_FILTER);
    }

    @Nullable
    public static MsJitterConfig create(final String sampleId, final String refGenomeFile, final RefGenomeVersion refGenVersion,
            final SequencingType sequencing, boolean usesDuplexUMIs, final String outputDir, final ConfigBuilder configBuilder)
    {
        if(!configBuilder.hasValue(JITTER_MSI_SITES_FILE))
            return null;

        return new MsJitterConfig(sampleId, refGenomeFile, refGenVersion, usesDuplexUMIs, outputDir, configBuilder);
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(JITTER_MSI_SITES_FILE, false, JITTER_MSI_SITES_FILE_DESC);
        configBuilder.addInteger(MIN_MAP_QUALITY, "Minimum mapping quality for an alignment to be used", DEFAULT_MIN_MAPPING_QUALITY);
        configBuilder.addInteger(JITTER_MAX_SITES_PER_TYPE, JITTER_MAX_SITES_PER_TYPE_DESC, DEFAULT_NUM_SITES_PER_TYPE);
        configBuilder.addDecimal(MAX_SINGLE_SITE_ALT_CONTRIBUTION, "Max percentage a single alt site can contribute", DEFAULT_MAX_SINGLE_SITE_ALT_CONTRIBUTION);
        configBuilder.addFlag(JITTER_WRITE_MSI_PLOTS, JITTER_WRITE_MSI_PLOTS_DESC);
        configBuilder.addFlag(JITTER_WRITE_SITE_FILE, JITTER_WRITE_SITE_FILE_DESC);

        configBuilder.addFlag(NO_BASE_QUAL_FILTER, "Include SBX medium base quals");
    }

    public static EnumSet<ConsensusType> consensusTypes(final MsJitterConfig config)
    {
        EnumSet<ConsensusType> consensusTypes = Sets.newEnumSet(Collections.emptyList(), ConsensusType.class);

        consensusTypes.add(ConsensusType.NONE);

        if(SEQUENCING_TYPE == ILLUMINA && config.UsesDuplexUMIs)
        {
            consensusTypes.add(ConsensusType.SINGLE);
            consensusTypes.add(ConsensusType.DUAL);
        }
        else if(SEQUENCING_TYPE == SBX)
        {
            consensusTypes.add(ConsensusType.SINGLE);
            consensusTypes.add(ConsensusType.DUAL);
        }
        else if(SEQUENCING_TYPE == ULTIMA)
        {
            consensusTypes.add(ConsensusType.DUAL);
        }
        else if(SEQUENCING_TYPE == BIOMODAL)
        {
            consensusTypes.add(ConsensusType.DUAL);
        }

        return consensusTypes;
    }

}
