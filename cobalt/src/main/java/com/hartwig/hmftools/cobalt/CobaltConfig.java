package com.hartwig.hmftools.cobalt;

import static com.hartwig.hmftools.cobalt.CobaltConstants.DEFAULT_MIN_MAPPING_QUALITY;
import static com.hartwig.hmftools.cobalt.CobaltConstants.DEFAULT_PCF_GAMMA;
import static com.hartwig.hmftools.common.genome.gc.GCProfileFactory.GC_PROFILE;
import static com.hartwig.hmftools.common.genome.gc.GCProfileFactory.GC_PROFILE_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME_CFG_DESC;
import static com.hartwig.hmftools.common.bam.BamUtils.addValidationStringencyOption;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE_BAM;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE_BAM_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_BAM;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_BAM_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkCreateOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;

import com.hartwig.hmftools.common.bam.BamUtils;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.ValidationStringency;

public class CobaltConfig
{
    public enum Mode
    {
        TUMOR_GERMLINE,
        TUMOR_ONLY,
        GERMLIHE_ONLY
    }

    private static final String TUMOR_ONLY_DIPLOID_BED = "tumor_only_diploid_bed";
    private static final String MIN_MAPPING_QUALITY = "min_quality";
    private static final String PCF_GAMMA = "pcf_gamma";
    private static final String TARGET_REGION_NORM_FILE = "target_region";
    private static final String INCLUDE_DUPLICATES = "include_duplicates";


    public final String ReferenceId;
    public final String ReferenceBamPath;

    public final String TumorId;
    public final String TumorBamPath;

    public final String RefGenomePath;
    public final String GcProfilePath;

    public final int Threads;
    public final String OutputDir;

    public final int MinMappingQuality;
    public final int PcfGamma;

    public final ValidationStringency BamStringency;
    public final boolean IncludeDuplicates;

    public final String TumorOnlyDiploidBed;
    public final String TargetRegionPath;

    public static final Logger CB_LOGGER = LogManager.getLogger(CobaltConfig.class);

    public CobaltConfig(final ConfigBuilder configBuilder)
    {
        TumorId = configBuilder.getValue(TUMOR);
        TumorBamPath = configBuilder.getValue(TUMOR_BAM);

        ReferenceId = configBuilder.getValue(REFERENCE);
        ReferenceBamPath = configBuilder.getValue(REFERENCE_BAM);

        GcProfilePath = configBuilder.getValue(GC_PROFILE);

        TumorOnlyDiploidBed = configBuilder.getValue(TUMOR_ONLY_DIPLOID_BED);
        TargetRegionPath = configBuilder.getValue(TARGET_REGION_NORM_FILE);
        RefGenomePath = configBuilder.getValue(REF_GENOME);
        
        MinMappingQuality = configBuilder.getInteger(MIN_MAPPING_QUALITY);
        PcfGamma = configBuilder.getInteger(PCF_GAMMA);
        IncludeDuplicates = configBuilder.hasFlag(INCLUDE_DUPLICATES);

        BamStringency = BamUtils.validationStringency(configBuilder);
        OutputDir = parseOutputDir(configBuilder);
        Threads = parseThreads(configBuilder);

    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(TUMOR, TUMOR_DESC);
        configBuilder.addPath(TUMOR_BAM, false, TUMOR_BAM_DESC);

        configBuilder.addConfigItem(REFERENCE, REFERENCE_DESC);
        configBuilder.addPath(REFERENCE_BAM, false, REFERENCE_BAM_DESC);

        configBuilder.addPath(GC_PROFILE, true, GC_PROFILE_DESC);
        configBuilder.addPath(REF_GENOME, false, REF_GENOME_CFG_DESC + ", required when using CRAM files");
        configBuilder.addPath(TUMOR_ONLY_DIPLOID_BED, false, "Diploid regions for tumor-only mode");
        configBuilder.addPath(TARGET_REGION_NORM_FILE, false, "Targeted regions normalisation file");

        configBuilder.addInteger(MIN_MAPPING_QUALITY, "Min map quality", DEFAULT_MIN_MAPPING_QUALITY);
        configBuilder.addInteger(PCF_GAMMA, "Gamma value for copy number PCF", DEFAULT_PCF_GAMMA);
        configBuilder.addFlag(INCLUDE_DUPLICATES, "Include duplicate reads in depth counts");

        addOutputDir(configBuilder);
        addThreadOptions(configBuilder);
        addValidationStringencyOption(configBuilder);
        addLoggingOptions(configBuilder);
    }

    public void validate() throws Exception
    {
        if(ReferenceId == null)
        {
            if(ReferenceBamPath != null)
            {
                throw new Exception(String.format("%s option not allowed in tumor only mode", REFERENCE_BAM));
            }
            if(TumorOnlyDiploidBed == null)
            {
                throw new Exception(String.format("missing required option %s in tumor only mode", TUMOR_ONLY_DIPLOID_BED));
            }
        }
        else if(TumorOnlyDiploidBed != null)
        {
            throw new Exception(String.format("%s option is only allowed in tumor only mode", TUMOR_ONLY_DIPLOID_BED));
        }

        if(TumorId == null)
        {
            if(TumorBamPath != null)
            {
                throw new Exception(String.format("%s option not allowed in germline only mode", TUMOR_BAM));
            }
        }

        if(!checkCreateOutputDir(OutputDir))
        {
            throw new Exception(String.format("failed to create output directory(%s)", OutputDir));
        }

        if(GcProfilePath.endsWith(".gz"))
        {
            throw new Exception(String.format("invalid GC-profile file(%s), must be uncompressed", GcProfilePath));
        }
    }

    public Mode mode()
    {
        if(ReferenceId != null && TumorId == null)
            return Mode.GERMLIHE_ONLY;

        if(ReferenceId == null && TumorId != null)
            return Mode.TUMOR_ONLY;

        return Mode.TUMOR_GERMLINE;
    }
}
