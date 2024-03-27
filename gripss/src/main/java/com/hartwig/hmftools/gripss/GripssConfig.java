package com.hartwig.hmftools.gripss;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.region.SpecificRegions.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.region.SpecificRegions.loadSpecificChromsomes;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.gripss.RepeatMaskAnnotations;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;
import com.hartwig.hmftools.gripss.common.SvData;
import com.hartwig.hmftools.gripss.filters.FilterConstants;
import com.hartwig.hmftools.gripss.filters.HotspotCache;
import com.hartwig.hmftools.gripss.filters.TargetRegions;
import com.hartwig.hmftools.gripss.pon.PonCache;

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

    private static final String VCF_FILE = "vcf";
    private static final String GERMLINE = "germline";

    public static final Logger GR_LOGGER = LogManager.getLogger(GripssApplication.class);

    public static final String APP_NAME = "Gripss";

    public GripssConfig(final ConfigBuilder configBuilder)
    {
        SampleId = configBuilder.getValue(SAMPLE);
        ReferenceId = configBuilder.getValue(REFERENCE, "");
        GermlineMode = configBuilder.hasFlag(GERMLINE);
        OutputDir = parseOutputDir(configBuilder);
        OutputId = configBuilder.getValue(OUTPUT_ID);

        VcfFile = configBuilder.getValue(VCF_FILE);
        RefGenVersion = RefGenomeVersion.from(configBuilder);

        RestrictedChromosomes = loadSpecificChromsomes(configBuilder);
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

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(SAMPLE, true, SAMPLE_DESC);
        configBuilder.addConfigItem(REFERENCE, REFERENCE_DESC);
        configBuilder.addPath(VCF_FILE, true, "Path to the GRIDSS structural variant VCF file");
        configBuilder.addFlag(GERMLINE, "Run in germline mode");

        addOutputOptions(configBuilder);
        ConfigUtils.addLoggingOptions(configBuilder);
        addRefGenomeConfig(configBuilder, true);

        addSpecificChromosomesRegionsConfig(configBuilder);

        PonCache.addConfig(configBuilder);
        HotspotCache.addConfig(configBuilder);
        FilterConstants.addConfig(configBuilder);
        RepeatMaskAnnotations.addConfig(configBuilder);
        TargetRegions.addConfig(configBuilder);
    }
}
