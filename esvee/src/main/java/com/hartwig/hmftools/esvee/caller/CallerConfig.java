package com.hartwig.hmftools.esvee.caller;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION_CFG_DESC;
import static com.hartwig.hmftools.common.region.SpecificRegions.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.region.SpecificRegions.loadSpecificChromsomes;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.caller.annotation.PonCache.ARTEFACT_PON_BED_SGL_FILE;
import static com.hartwig.hmftools.esvee.caller.annotation.PonCache.ARTEFACT_PON_BED_SV_FILE;
import static com.hartwig.hmftools.esvee.common.FileCommon.DEPTH_VCF_SUFFIX;
import static com.hartwig.hmftools.esvee.common.FileCommon.INPUT_VCF;
import static com.hartwig.hmftools.esvee.common.FileCommon.INPUT_VCF_DESC;
import static com.hartwig.hmftools.esvee.common.FileCommon.PREP_DIR;
import static com.hartwig.hmftools.esvee.common.FileCommon.PREP_DIR_DESC;
import static com.hartwig.hmftools.esvee.common.FileCommon.formEsveeInputFilename;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.sv.RepeatMaskAnnotations;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.esvee.caller.annotation.PonCache;

public class CallerConfig
{
    // run config
    public final String TumorId;
    public final String ReferenceId;
    public final RefGenomeVersion RefGenVersion;

    public final String VcfFile;
    public final String PrepDir;

    public final String OutputDir;
    public final String OutputId;
    public final List<String> SpecificChromosomes;

    public final int ManualRefDepth;
    public final boolean WriteBreakendTsv;

    public static final String MANUAL_REF_DEPTH = "manual_ref_depth";
    public static final String WRITE_BREAKEND_TSV = "write_breakend_tsv";

    public CallerConfig(final ConfigBuilder configBuilder)
    {
        TumorId = configBuilder.getValue(TUMOR);
        ReferenceId = configBuilder.getValue(REFERENCE);

        OutputDir = parseOutputDir(configBuilder);
        OutputId = configBuilder.getValue(OUTPUT_ID);

        if(configBuilder.hasValue(INPUT_VCF))
        {
            VcfFile = configBuilder.getValue(INPUT_VCF);
        }
        else
        {
            String fileSampleId = TumorId != null ? TumorId : ReferenceId;
            VcfFile = formEsveeInputFilename(OutputDir, fileSampleId, DEPTH_VCF_SUFFIX, OutputId);
        }

        PrepDir = configBuilder.hasValue(PREP_DIR) ? checkAddDirSeparator(configBuilder.getValue(PREP_DIR)) : OutputDir;

        RefGenVersion = RefGenomeVersion.from(configBuilder);

        SpecificChromosomes = loadSpecificChromsomes(configBuilder);

        ManualRefDepth = configBuilder.getInteger(MANUAL_REF_DEPTH);
        WriteBreakendTsv = configBuilder.hasFlag(WRITE_BREAKEND_TSV);
    }

    public boolean hasTumor() { return TumorId != null; }
    public boolean hasReference() { return ReferenceId != null; }
    public boolean germlineOnly() { return ReferenceId != null && TumorId == null;}

    public String fileSampleId() { return TumorId != null ? TumorId : ReferenceId; }

    public boolean isValid()
    {
        if(TumorId == null && ReferenceId == null)
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
        if(!SpecificChromosomes.isEmpty() && !SpecificChromosomes.contains(sv.chromosomeStart())
        && (sv.type() == SGL || !SpecificChromosomes.contains(sv.chromosomeEnd())))
        {
            return true;
        }

        return false;
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        if(!configBuilder.isRegistered(TUMOR))
            configBuilder.addConfigItem(TUMOR, TUMOR_DESC);

        if(!configBuilder.isRegistered(REFERENCE))
            configBuilder.addConfigItem(REFERENCE, REFERENCE_DESC);

        configBuilder.addPath(INPUT_VCF, false, INPUT_VCF_DESC);
        configBuilder.addPaths(PREP_DIR, false, PREP_DIR_DESC);
        configBuilder.addInteger(MANUAL_REF_DEPTH, "Manually set ref depth for testing", 0);

        addOutputOptions(configBuilder);
        addLoggingOptions(configBuilder);
        configBuilder.addConfigItem(REF_GENOME_VERSION, true, REF_GENOME_VERSION_CFG_DESC);

        addSpecificChromosomesRegionsConfig(configBuilder);

        PonCache.addConfig(configBuilder);

        configBuilder.addPath(ARTEFACT_PON_BED_SV_FILE, false, "Additional artefact SV PON file");
        configBuilder.addPath(ARTEFACT_PON_BED_SGL_FILE, false, "Additional artefact SGL PON file");

        configBuilder.addFlag(WRITE_BREAKEND_TSV, "Rewrite a breakend TSV for with additional caller annotations");

        HotspotCache.addConfig(configBuilder);
        FilterConstants.addConfig(configBuilder);
        RepeatMaskAnnotations.addConfig(configBuilder);
        TargetRegions.addConfig(configBuilder);
    }
}
