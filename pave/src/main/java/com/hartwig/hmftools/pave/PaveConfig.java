package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.region.SpecificRegions.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.region.SpecificRegions.loadSpecificRegions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelConfig;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.pave.annotation.Blacklistings;
import com.hartwig.hmftools.pave.annotation.ClinvarAnnotation;
import com.hartwig.hmftools.pave.annotation.GnomadAnnotation;
import com.hartwig.hmftools.pave.annotation.Mappability;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class PaveConfig
{
    public final String SampleId;

    public final String VcfFile;
    public final RefGenomeVersion RefGenVersion;

    public final String OutputDir;
    public final String OutputVcfFile;

    public final boolean WriteTranscriptFile;
    public final boolean OnlyCanonical;
    public final boolean ReadPassOnly;
    public final boolean WritePassOnly;


    // PON related threholds
    // public final double mPonFilterThreshold;

    public final boolean SetReportable;
    public final List<ChrBaseRegion> SpecificRegions;
    public final int Threads;

    public static final String VCF_FILE = "vcf_file";
    private static final String OUTPUT_VCF_FILE = "output_vcf_file";

    public static final String PON_FILE = "pon_file";
    public static final String PON_FILTERS = "pon_filters";
    public static final String PON_ARTEFACTS_FILE = "pon_artefact_file";

    // optional and debugging config
    private static final String ONLY_CANONCIAL = "only_canonical";
    private static final String READ_PASS_ONLY = "read_pass_only";
    private static final String WRITE_PASS_ONLY = "write_pass_only";
    private static final String WRITE_TRANSCRIPT_DATA = "write_transcript_data";
    private static final String SET_REPORTABLE = "set_reportable";

    public static final Logger PV_LOGGER = LogManager.getLogger(PaveConfig.class);

    public PaveConfig(final ConfigBuilder configBuilder)
    {
        SampleId = configBuilder.getValue(SAMPLE);
        VcfFile = configBuilder.getValue(VCF_FILE);

        RefGenVersion = RefGenomeVersion.from(configBuilder);

        OutputVcfFile = configBuilder.getValue(OUTPUT_VCF_FILE);

        WriteTranscriptFile = SampleId != null && configBuilder.hasFlag(WRITE_TRANSCRIPT_DATA);
        OnlyCanonical = configBuilder.hasFlag(ONLY_CANONCIAL);
        ReadPassOnly = configBuilder.hasFlag(READ_PASS_ONLY);
        WritePassOnly = configBuilder.hasFlag(WRITE_PASS_ONLY);
        SetReportable = configBuilder.hasFlag(SET_REPORTABLE);
        Threads = parseThreads(configBuilder);

        SpecificRegions = Lists.newArrayList();

        try
        {
            SpecificRegions.addAll(loadSpecificRegions(configBuilder));
        }
        catch(Exception e)
        {
            PV_LOGGER.error("failed to load specific regions");
        }

        if(configBuilder.hasValue(OUTPUT_DIR))
        {
            OutputDir = parseOutputDir(configBuilder);
        }
        else
        {
            String vcfFile = OutputVcfFile != null ? OutputVcfFile : VcfFile;
            Path vcfDir = Paths.get(vcfFile).getParent();
            OutputDir = vcfDir != null ? checkAddDirSeparator(vcfDir.toString()) : "./";
        }
    }

    public boolean isValid()
    {
        if(VcfFile == null)
            return false;

        return true;
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(SAMPLE, true, SAMPLE_DESC);
        configBuilder.addPath(VCF_FILE, true, "VCF input file");
        configBuilder.addConfigItem(
                OUTPUT_VCF_FILE, false, "Option VCF output file, otherwise will append 'pave' suffix to input filename");

        addRefGenomeConfig(configBuilder, true);
        addEnsemblDir(configBuilder, true);
        DriverGenePanelConfig.addGenePanelOption(configBuilder, false);
        configBuilder.addPath(PON_FILE, false, "PON entries");
        configBuilder.addPath(PON_ARTEFACTS_FILE, false, "PON artefacts to filter");
        configBuilder.addConfigItem(PON_FILTERS, "PON filters per tier, format: TIER:MAX_SAMPLES:MAX_COUNT separated by ';'");

        configBuilder.addFlag(WRITE_TRANSCRIPT_DATA, "Write variant impacts per transcript to TSV");
        configBuilder.addFlag(ONLY_CANONCIAL, "Only check canonical transcripts");
        configBuilder.addFlag(READ_PASS_ONLY, "Filter incoming variants to PASS only");
        configBuilder.addFlag(WRITE_PASS_ONLY, "Only annotate passing variants");
        configBuilder.addFlag(SET_REPORTABLE, "Set reportable and hotspot flags");

        GnomadAnnotation.addConfig(configBuilder);
        Mappability.addConfig(configBuilder);
        ClinvarAnnotation.addConfig(configBuilder);
        Blacklistings.addConfig(configBuilder);
        addThreadOptions(configBuilder);
        addSpecificChromosomesRegionsConfig(configBuilder);

        addOutputDir(configBuilder);
        ConfigUtils.addLoggingOptions(configBuilder);
    }
}
