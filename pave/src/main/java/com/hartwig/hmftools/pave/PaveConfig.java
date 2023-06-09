package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION_CFG_DESC;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.loadSpecificRegions;

import java.nio.file.Paths;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelConfig;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.pave.annotation.Blacklistings;
import com.hartwig.hmftools.pave.annotation.ClinvarAnnotation;
import com.hartwig.hmftools.pave.annotation.GnomadAnnotation;
import com.hartwig.hmftools.pave.annotation.Mappability;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class PaveConfig
{
    public final String SampleId;

    public final String VcfFile;
    public final RefGenomeVersion RefGenVersion;

    public final String OutputDir;
    public final String OutputVcfFile;

    public final boolean WriteTranscriptFile;
    public final boolean WriteDiffs;
    public final boolean OnlyCanonical;
    public final boolean ReadPassOnly;
    public final boolean WritePassOnly;
    public final boolean SetReportable;
    public final List<ChrBaseRegion> SpecificRegions;

    private static final String SAMPLE = "sample";
    private static final String VCF_FILE = "vcf_file";
    private static final String OUTPUT_VCF_FILE = "output_vcf_file";

    public static final String PON_FILE = "pon_file";
    public static final String PON_FILTERS = "pon_filters";
    public static final String PON_ARTEFACTS_FILE = "pon_artefact_file";

    // optional and debugging config
    private static final String ONLY_CANONCIAL = "only_canonical";
    private static final String READ_PASS_ONLY = "read_pass_only";
    private static final String WRITE_PASS_ONLY = "write_pass_only";
    private static final String WRITE_DIFFS = "write_diffs";
    private static final String WRITE_TRANSCRIPT_DATA = "write_transcript_data";
    private static final String SET_REPORTABLE = "set_reportable";

    public static final Logger PV_LOGGER = LogManager.getLogger(PaveConfig.class);

    public PaveConfig(final CommandLine cmd)
    {
        SampleId = cmd.getOptionValue(SAMPLE);
        VcfFile = cmd.getOptionValue(VCF_FILE);

        RefGenVersion = RefGenomeVersion.from(cmd);

        OutputVcfFile = cmd.getOptionValue(OUTPUT_VCF_FILE);

        WriteTranscriptFile = SampleId != null && cmd.hasOption(WRITE_TRANSCRIPT_DATA);
        WriteDiffs = cmd.hasOption(WRITE_DIFFS);
        OnlyCanonical = cmd.hasOption(ONLY_CANONCIAL);
        ReadPassOnly = cmd.hasOption(READ_PASS_ONLY);
        WritePassOnly = cmd.hasOption(WRITE_PASS_ONLY);
        SetReportable = cmd.hasOption(SET_REPORTABLE);

        SpecificRegions = Lists.newArrayList();

        try
        {
            SpecificRegions.addAll(loadSpecificRegions(cmd));
        }
        catch(Exception e)
        {
            PV_LOGGER.error("failed to load specific regions");
        }

        if(cmd.hasOption(OUTPUT_DIR))
        {
            OutputDir = parseOutputDir(cmd);
        }
        else
        {
            OutputDir = checkAddDirSeparator(Paths.get(VcfFile).getParent().toString());
        }
    }

    public boolean isValid()
    {
        if(VcfFile == null)
            return false;

        return true;
    }

    @NotNull
    public static Options createOptions()
    {
        Options options = new Options();
        options.addOption(SAMPLE, true, "Name of sample");
        options.addOption(VCF_FILE, true, "VCF input file");
        options.addOption(OUTPUT_VCF_FILE, true, "Option VCF output file, otherwise will append 'pave' suffix to input filename");

        options.addOption(REF_GENOME, true, REF_GENOME_CFG_DESC);
        options.addOption(REF_GENOME_VERSION, true, REF_GENOME_VERSION_CFG_DESC);
        addEnsemblDir(options);
        DriverGenePanelConfig.addGenePanelOption(false, options);
        options.addOption(PON_FILE, true, "PON entries");
        options.addOption(PON_ARTEFACTS_FILE, true, "PON artefacts to filter");
        options.addOption(PON_FILTERS, true, "PON filters per tier, format: TIER:MAX_SAMPLES:MAX_COUNT separated by ';'");

        options.addOption(WRITE_DIFFS, false, "Only write transcript diffs to CSV file");
        options.addOption(WRITE_TRANSCRIPT_DATA, false, "Write variant impacts per transcript to TSV");
        options.addOption(ONLY_CANONCIAL, false, "Only check canonical transcripts");
        options.addOption(READ_PASS_ONLY, false, "Filter incoming variants to PASS only");
        options.addOption(WRITE_PASS_ONLY, false, "Only annotate passing variants");
        options.addOption(SET_REPORTABLE, false, "Set reportable and hotspot flags");

        GnomadAnnotation.addCmdLineArgs(options);
        Mappability.addCmdLineArgs(options);
        ClinvarAnnotation.addCmdLineArgs(options);
        Blacklistings.addCmdLineArgs(options);
        addSpecificChromosomesRegionsConfig(options);

        addOutputDir(options);
        addLoggingOptions(options);

        return options;
    }
}
