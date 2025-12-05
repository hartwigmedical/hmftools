package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.region.SpecificRegions.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.region.SpecificRegions.loadSpecificRegions;
import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.perf.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigItem.enumValueSelectionAsStr;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.VCF_ZIP_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.pathFromFile;
import static com.hartwig.hmftools.common.variant.pon.PonCache.PON_FILE;
import static com.hartwig.hmftools.common.variant.pon.PonCache.PON_FILTERS;
import static com.hartwig.hmftools.pave.FilterType.ALL;
import static com.hartwig.hmftools.pave.FilterType.PASS;
import static com.hartwig.hmftools.pave.annotation.GnomadAnnotation.GNOMAD_NO_FILTER;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.driver.panel.DriverGenePanelConfig;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.file.FileReaderUtils;
import com.hartwig.hmftools.common.variant.pon.GnomadCache;
import com.hartwig.hmftools.pave.annotation.Blacklistings;
import com.hartwig.hmftools.pave.annotation.ClinvarAnnotation;
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
    public final FilterType Filter;
    public final boolean WritePassOnly;
    public final boolean WriteDetailed;

    public final boolean SetReportable;
    public final List<ChrBaseRegion> SpecificRegions;
    public final int Threads;

    private static final String INPUT_VCF = "input_vcf";
    private static final String OUTPUT_VCF = "output_vcf";

    public static final String PON_ARTEFACTS_FILE = "pon_artefact_file";

    // optional and debugging config
    private static final String ONLY_CANONCIAL = "only_canonical";
    private static final String FILTER_TYPE = "filter_type";
    private static final String PROCESS_NON_PASS = "process_non_pass"; // to be deprecated
    private static final String WRITE_PASS_ONLY = "write_pass_only";
    private static final String WRITE_TRANSCRIPT_DATA = "write_transcript_data";
    private static final String WRITE_DETAILED = "write_detailed";
    private static final String SET_REPORTABLE = "set_reportable";

    public static final Logger PV_LOGGER = LogManager.getLogger(PaveConfig.class);

    public PaveConfig(final ConfigBuilder configBuilder)
    {
        SampleId = configBuilder.getValue(SAMPLE);

        VcfFile = configBuilder.getValue(INPUT_VCF);
        OutputVcfFile = configBuilder.getValue(OUTPUT_VCF);

        RefGenVersion = RefGenomeVersion.from(configBuilder);

        WriteTranscriptFile = SampleId != null && configBuilder.hasFlag(WRITE_TRANSCRIPT_DATA);
        OnlyCanonical = configBuilder.hasFlag(ONLY_CANONCIAL);

        if(configBuilder.hasValue(FILTER_TYPE))
        {
            Filter = FilterType.valueOf(configBuilder.getValue(FILTER_TYPE));
        }
        else if(configBuilder.hasValue(PROCESS_NON_PASS))
        {
            Filter = ALL;
        }
        else
        {
            Filter = PASS;
        }

        WritePassOnly = configBuilder.hasFlag(WRITE_PASS_ONLY);
        WriteDetailed = configBuilder.hasFlag(WRITE_DETAILED);
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
            OutputDir = pathFromFile(vcfFile);
        }
    }

    public boolean requireIndex() { return VcfFile.endsWith(VCF_ZIP_EXTENSION); }

    public boolean isValid()
    {
        if(VcfFile == null)
            return false;

        return true;
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(SAMPLE, true, SAMPLE_DESC);

        configBuilder.addPath(INPUT_VCF, true, "VCF input file");
        configBuilder.addConfigItem(
                OUTPUT_VCF, false, "VCF output file (optional), otherwise will append 'pave' suffix to input filename");

        addRefGenomeConfig(configBuilder, true);
        addEnsemblDir(configBuilder, true);
        DriverGenePanelConfig.addGenePanelOption(configBuilder, false);

        configBuilder.addPath(PON_FILE, false, "PON entries");
        configBuilder.addPath(PON_ARTEFACTS_FILE, false, "PON artefacts to filter");
        configBuilder.addConfigItem(PON_FILTERS, "PON filters per tier, format: TIER:MAX_SAMPLES:MAX_COUNT separated by ';'");

        configBuilder.addFlag(WRITE_TRANSCRIPT_DATA, "Write variant impacts per transcript to TSV");
        configBuilder.addFlag(ONLY_CANONCIAL, "Only check canonical transcripts");

        configBuilder.addFlag(PROCESS_NON_PASS, "Process all variants from Sage");
        configBuilder.addConfigItem(FILTER_TYPE, false, enumValueSelectionAsStr(FilterType.values(), "Variant filters"));

        configBuilder.addFlag(WRITE_PASS_ONLY, "Only annotate passing variants");
        configBuilder.addFlag(SET_REPORTABLE, "Set reportable and hotspot flags");
        configBuilder.addFlag(WRITE_DETAILED, "Write detailed transcript impact info");

        GnomadCache.addConfig(configBuilder);
        configBuilder.addFlag(GNOMAD_NO_FILTER, "No Gnomad filter is applied");
        Mappability.addConfig(configBuilder);
        ClinvarAnnotation.addConfig(configBuilder);
        Blacklistings.addConfig(configBuilder);
        addThreadOptions(configBuilder);
        addSpecificChromosomesRegionsConfig(configBuilder);

        addOutputDir(configBuilder);
        ConfigUtils.addLoggingOptions(configBuilder);
    }
}
