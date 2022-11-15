package com.hartwig.hmftools.gripss;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.SPECIFIC_CHROMOSOMES;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.SPECIFIC_CHROMOSOMES_DESC;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.loadSpecificChromsomes;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.gripss.common.SvData;
import com.hartwig.hmftools.gripss.filters.FilterConstants;
import com.hartwig.hmftools.gripss.filters.HotspotCache;
import com.hartwig.hmftools.gripss.filters.TargetRegions;
import com.hartwig.hmftools.gripss.rm.RepeatMaskAnnotations;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class GripssConfig
{
    // run config
    public final String SampleId;
    public final String ReferenceId;
    public final RefGenomeVersion RefGenVersion;
    public final String VcfFile;
    public final boolean GermlineMode;

    public final String OutputDir;
    public final String OutputId;
    public final List<String> RestrictedChromosomes;

    public static final String SAMPLE = "sample";
    private static final String REFERENCE = "reference";
    private static final String VCF_FILE = "vcf";
    private static final String GERMLINE = "germline";

    public static final Logger GR_LOGGER = LogManager.getLogger(GripssApplication.class);

    public GripssConfig(final CommandLine cmd)
    {
        SampleId = cmd.getOptionValue(SAMPLE, "");
        ReferenceId = cmd.getOptionValue(REFERENCE, "");
        GermlineMode = cmd.hasOption(GERMLINE);
        OutputDir = parseOutputDir(cmd);
        OutputId = cmd.getOptionValue(OUTPUT_ID);

        VcfFile = cmd.getOptionValue(VCF_FILE, "");
        RefGenVersion = RefGenomeVersion.from(cmd.getOptionValue(REF_GENOME_VERSION, V37.toString()));

        RestrictedChromosomes = loadSpecificChromsomes(cmd);
    }

    public GripssConfig(
            final String sampleId, final String referenceId, final RefGenomeVersion refGenVersion, final String vcfFile)
    {
        SampleId = sampleId;
        ReferenceId = referenceId;
        RefGenVersion = refGenVersion;
        GermlineMode = false;
        VcfFile = vcfFile;
        OutputDir = null;
        OutputId = null;
        RestrictedChromosomes = Lists.newArrayList();
    }

    public boolean tumorOnly() { return ReferenceId.isEmpty(); }

    public boolean isValid()
    {
        if(SampleId.isEmpty())
        {
            GR_LOGGER.error("missing sample config");
            return false;
        }

        if(VcfFile.isEmpty() || !Files.exists(Paths.get(VcfFile)))
        {
            GR_LOGGER.error("missing or invalid VCF file({})", VcfFile);
            return false;
        }

        return true;
    }

    public boolean excludeVariant(final SvData sv)
    {
        // optionally filter out all but specified chromosomes
        if(!RestrictedChromosomes.isEmpty() && !RestrictedChromosomes.contains(sv.chromosomeStart())
        && (sv.type() == SGL || !RestrictedChromosomes.contains(sv.chromosomeEnd())))
        {
            return true;
        }

        return false;
    }

    public static void addCommandLineOptions(final Options options)
    {
        options.addOption(SAMPLE, true, "Name of the tumor sample");
        options.addOption(REFERENCE, true, "Optional, name of the reference sample");
        options.addOption(VCF_FILE, true, "Path to the GRIDSS structural variant VCF file");
        options.addOption(GERMLINE, false, "Run in germline mode");
        options.addOption(REF_GENOME, true, REF_GENOME_CFG_DESC);
        options.addOption(REF_GENOME_VERSION, true, REF_GENOME_VERSION_CFG_DESC);

        addOutputOptions(options);
        addLoggingOptions(options);
        addRefGenomeConfig(options);

        options.addOption(SPECIFIC_CHROMOSOMES, true, SPECIFIC_CHROMOSOMES_DESC);

        PonCache.addCmdLineArgs(options);
        HotspotCache.addCmdLineArgs(options);
        FilterConstants.addCmdLineArgs(options);
        RepeatMaskAnnotations.addCmdLineArgs(options);
        TargetRegions.addCmdLineArgs(options);
    }
}
