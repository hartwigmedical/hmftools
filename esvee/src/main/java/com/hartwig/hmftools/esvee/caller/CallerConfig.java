package com.hartwig.hmftools.esvee.caller;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION_CFG_DESC;
import static com.hartwig.hmftools.common.region.SpecificRegions.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.region.SpecificRegions.loadSpecificChromsomes;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.common.FileCommon.DEPTH_VCF_SUFFIX;
import static com.hartwig.hmftools.esvee.common.FileCommon.INPUT_VCF;
import static com.hartwig.hmftools.esvee.common.FileCommon.INPUT_VCF_DESC;
import static com.hartwig.hmftools.esvee.common.FileCommon.formEsveeInputFilename;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.gripss.RepeatMaskAnnotations;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.esvee.assembly.output.WriteType;
import com.hartwig.hmftools.esvee.caller.annotation.PonCache;

public class CallerConfig
{
    // run config
    public final String SampleId;
    public final String ReferenceId;
    public final boolean GermlineOnly;
    public final RefGenomeVersion RefGenVersion;
    public final String VcfFile;

    public final String OutputDir;
    public final String OutputId;
    public final List<String> RestrictedChromosomes;

    public CallerConfig(final ConfigBuilder configBuilder)
    {
        SampleId = configBuilder.getValue(SAMPLE);
        ReferenceId = configBuilder.getValue(REFERENCE);

        GermlineOnly = ReferenceId != null && SampleId == null;

        OutputDir = parseOutputDir(configBuilder);
        OutputId = configBuilder.getValue(OUTPUT_ID);

        if(configBuilder.hasValue(INPUT_VCF))
            VcfFile = configBuilder.getValue(INPUT_VCF);
        else
            VcfFile = formEsveeInputFilename(OutputDir, SampleId, DEPTH_VCF_SUFFIX, OutputId);

        RefGenVersion = RefGenomeVersion.from(configBuilder);

        RestrictedChromosomes = loadSpecificChromsomes(configBuilder);
    }

    public boolean hasTumor() { return SampleId != null; }
    public boolean hasReference() { return ReferenceId != null; }

    public boolean isValid()
    {
        if(SampleId.isEmpty())
        {
            SV_LOGGER.error("missing sample config");
            return false;
        }

        if(VcfFile.isEmpty() || !Files.exists(Paths.get(VcfFile)))
        {
            SV_LOGGER.error("missing or invalid VCF file({})", VcfFile);
            return false;
        }

        return true;
    }

    public boolean excludeVariant(final Variant sv)
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
        configBuilder.addPath(INPUT_VCF, false, INPUT_VCF_DESC);

        addOutputOptions(configBuilder);
        addLoggingOptions(configBuilder);
        configBuilder.addConfigItem(REF_GENOME_VERSION, true, REF_GENOME_VERSION_CFG_DESC);

        addSpecificChromosomesRegionsConfig(configBuilder);

        PonCache.addConfig(configBuilder);
        HotspotCache.addConfig(configBuilder);
        FilterConstants.addConfig(configBuilder);
        RepeatMaskAnnotations.addConfig(configBuilder);
        TargetRegions.addConfig(configBuilder);
    }
}
