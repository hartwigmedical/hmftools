package com.hartwig.hmftools.cobalt;

import static com.hartwig.hmftools.cobalt.CobaltConstants.DEFAULT_MIN_MAPPING_QUALITY;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME_CFG_DESC;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkCreateOutputDir;

import java.io.File;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParameterException;

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

    private static final int DEFAULT_THREADS = 4;

    public static final String TUMOR = "-tumor";
    public static final String REFERENCE = "-reference";

    private static final String TUMOR_ONLY_DIPLOID_BED = "-tumor_only_diploid_bed";
    private static final String REFERENCE_BAM = "-reference_bam";
    private static final String TUMOR_BAM = "-tumor_bam";
    private static final String GC_PROFILE = "-gc_profile";
    private static final String MIN_MAPPING_QUALITY = "-min_quality";
    private static final String VALIDATION_STRINGENCY = "-validation_stringency";

    @Parameter(names = "-threads", description = "Number of threads")
    public int ThreadCount = DEFAULT_THREADS;

    @Parameter(names = MIN_MAPPING_QUALITY, description = "Min quality")
    public int MinMappingQuality = DEFAULT_MIN_MAPPING_QUALITY;

    @Parameter(names = REFERENCE, description = "Name of reference sample")
    public String ReferenceId;

    @Parameter(names = REFERENCE_BAM, description = "Path to reference bam file")
    public String ReferenceBamPath;

    @Parameter(names = TUMOR, description = "Name of tumor sample")
    public String TumorId;

    @Parameter(names = TUMOR_BAM, description = "Path to tumor bam file")
    public String TumorBamPath;

    @Parameter(names = "-" + REF_GENOME, description = REF_GENOME_CFG_DESC)
    public String RefGenomePath;

    @Parameter(names = GC_PROFILE,
               required = true,
               description = "Location of GC Profile")
    public String GcProfilePath;

    @Parameter(names = "-" + OUTPUT_DIR,
               required = true,
               description = "Path to the output directory. This directory will be created if it does not already exist.")
    public String OutputDir;

    @Parameter(names = VALIDATION_STRINGENCY,
               description = "SAM validation strategy")
    public ValidationStringency Stringency = ValidationStringency.DEFAULT_STRINGENCY;

    @Parameter(names = TUMOR_ONLY_DIPLOID_BED,
               description = "Diploid regions for tumor-only mode")
    public String TumorOnlyDiploidBed;

    @Parameter(names = "-target_region",
               description = "Input tsv file for genome regions enrichment in targeted mode")
    public String TargetRegionPath;

    @Parameter(names = "-pcf_gamma",
               description = "Gamma value we pass to copynumber PCF")
    public int PcfGamma = CobaltConstants.DEFAULT_PCF_GAMMA;

    public static final Logger CB_LOGGER = LogManager.getLogger(CobaltConfig.class);

    public CobaltConfig()
    {
    }

    public void validate() throws ParameterException
    {
        if (ReferenceId == null)
        {
            if (ReferenceBamPath != null)
            {
                throw new ParameterException(String.format("%s option not allowed in tumor only mode", REFERENCE_BAM));
            }
            if (TumorOnlyDiploidBed == null)
            {
                throw new ParameterException(String.format("missing required option %s in tumor only mode", TUMOR_ONLY_DIPLOID_BED));
            }
            if (!new File(TumorOnlyDiploidBed).exists())
            {
                throw new ParameterException(String.format("diploid bed file %s does not exist", TumorOnlyDiploidBed));
            }
        }
        else if (TumorOnlyDiploidBed != null)
        {
            throw new ParameterException(String.format("%s option is only allowed in tumor only mode", TUMOR_ONLY_DIPLOID_BED));
        }

        if (TumorId == null)
        {
            if (TumorBamPath != null)
            {
                throw new ParameterException(String.format("%s option not allowed in germline only mode", TUMOR_BAM));
            }
        }

        if (!checkCreateOutputDir(OutputDir))
        {
            throw new ParameterException(String.format("failed to create output directory(%s)", OutputDir));
        }

        if (GcProfilePath.endsWith(".gz"))
        {
            throw new ParameterException(String.format("invalid GC-profile file(%s), must be uncompressed", GcProfilePath));
        }

        if (RefGenomePath != null && !(new File(RefGenomePath).exists()))
        {
            throw new ParameterException("Unable to locate ref genome file: " + RefGenomePath);
        }
    }

    public Mode mode()
    {
        if (ReferenceId != null && TumorId == null)
        {
            return Mode.GERMLIHE_ONLY;
        }
        if (ReferenceId == null && TumorId != null)
        {
            return Mode.TUMOR_ONLY;
        }
        return Mode.TUMOR_GERMLINE;
    }
}
