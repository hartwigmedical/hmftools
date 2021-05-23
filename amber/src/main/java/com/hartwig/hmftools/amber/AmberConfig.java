package com.hartwig.hmftools.amber;

import static com.hartwig.hmftools.amber.AmberConstants.DEFAULT_MAX_DEPTH_PERCENTAGE;
import static com.hartwig.hmftools.amber.AmberConstants.DEFAULT_MAX_HET_AF_PERCENTAGE;
import static com.hartwig.hmftools.amber.AmberConstants.DEFAULT_MIN_BASE_QUALITY;
import static com.hartwig.hmftools.amber.AmberConstants.DEFAULT_MIN_DEPTH_PERCENTAGE;
import static com.hartwig.hmftools.amber.AmberConstants.DEFAULT_MIN_HET_AF_PERCENTAGE;
import static com.hartwig.hmftools.amber.AmberConstants.DEFAULT_MIN_MAPPING_QUALITY;
import static com.hartwig.hmftools.amber.AmberConstants.DEFAULT_MIN_PARTITION;
import static com.hartwig.hmftools.amber.AmberConstants.DEFAULT_TUMOR_ONLY_MIN_SUPPORT;
import static com.hartwig.hmftools.amber.AmberConstants.DEFAULT_TUMOR_ONLY_MIN_VAF;
import static com.hartwig.hmftools.amber.AmberConstants.DEFAULT_TYPICAL_READ_DEPTH;
import static com.hartwig.hmftools.common.utils.ConfigUtils.getConfigValue;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkCreateOutputDir;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;

import static htsjdk.samtools.ValidationStringency.DEFAULT_STRINGENCY;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.ValidationStringency;

public class AmberConfig
{
    public final String TumorId;

    public final String BafLociPath;
    public final String TumorBamPath;
    public final List<String> ReferenceIds;
    public final List<String> ReferenceBamPath;
    public final String RefGenomePath;
    public final String OutputDir;

    public final boolean TumorOnly;
    public final int TumorOnlyMinSupport;
    public final double TumorOnlyMinVaf;

    public final int MinBaseQuality;
    public final int MinMappingQuality;
    public final double MinDepthPercent;
    public final double MaxDepthPercent;
    public final double MinHetAfPercent;
    public final double MaxHetAfPercent;
    public final ValidationStringency Stringency;

    public final int ThreadCount;

    public static final int DEFAULT_THREADS = 1;

    public static final String TUMOR = "tumor";
    public static final String BAF_LOCI = "loci";
    public static final String THREADS = "threads";
    public static final String REFERENCE = "reference";
    public static final String TUMOR_BAM = "tumor_bam";
    public static final String REF_GENOME = "ref_genome";
    public static final String OUTPUT_DIR = "output_dir";
    public static final String REFERENCE_BAM = "reference_bam";
    public static final String MIN_BASE_QUALITY = "min_base_quality";
    public static final String MIN_MAPPING_QUALITY = "min_mapping_quality";
    public static final String MIN_DEPTH_PERCENTAGE = "min_depth_percent";
    public static final String MAX_DEPTH_PERCENTAGE = "max_depth_percent";
    public static final String MIN_HET_AF_PERCENTAGE = "min_het_af_percent";
    public static final String MAX_HET_AF_PERCENTAGE = "max_het_af_percent";
    public static final String VALIDATION_STRINGENCY = "validation_stringency";

    public static final String TUMOR_ONLY = "tumor_only";
    public static final String TUMOR_ONLY_MIN_VAF = "tumor_only_min_vaf";
    public static final String TUMOR_ONLY_MIN_SUPPORT = "tumor_only_min_support";

    public static final Logger AMB_LOGGER = LogManager.getLogger(AmberConfig.class);

    public AmberConfig(final CommandLine cmd)
    {
        TumorId = cmd.getOptionValue(TUMOR);

        TumorOnly = cmd.hasOption(TUMOR_ONLY);

        ThreadCount = getConfigValue(cmd, THREADS, DEFAULT_THREADS);
        MinBaseQuality = getConfigValue(cmd, MIN_BASE_QUALITY, DEFAULT_MIN_BASE_QUALITY);
        MinMappingQuality = getConfigValue(cmd, MIN_MAPPING_QUALITY, DEFAULT_MIN_MAPPING_QUALITY);
        TumorOnlyMinSupport = getConfigValue(cmd, TUMOR_ONLY_MIN_SUPPORT, DEFAULT_TUMOR_ONLY_MIN_SUPPORT);

        MinDepthPercent = getConfigValue(cmd, MIN_DEPTH_PERCENTAGE, DEFAULT_MIN_DEPTH_PERCENTAGE);
        MaxDepthPercent = getConfigValue(cmd, MAX_DEPTH_PERCENTAGE, DEFAULT_MAX_DEPTH_PERCENTAGE);
        MinHetAfPercent = getConfigValue(cmd, MIN_HET_AF_PERCENTAGE, DEFAULT_MIN_HET_AF_PERCENTAGE);
        MaxHetAfPercent = getConfigValue(cmd, MAX_HET_AF_PERCENTAGE, DEFAULT_MAX_HET_AF_PERCENTAGE);
        TumorOnlyMinVaf = getConfigValue(cmd, TUMOR_ONLY_MIN_VAF, DEFAULT_TUMOR_ONLY_MIN_VAF);

        ReferenceIds = TumorOnly ? Collections.emptyList() : Arrays.asList(cmd.getOptionValue(REFERENCE).split(","));

        ReferenceBamPath = TumorOnly ? Collections.emptyList() : Arrays.asList(cmd.getOptionValue(REFERENCE_BAM).split(","));
        BafLociPath = cmd.getOptionValue(BAF_LOCI);
        TumorBamPath = cmd.getOptionValue(TUMOR_BAM);
        OutputDir = parseOutputDir(cmd);

        RefGenomePath = cmd.getOptionValue(REF_GENOME, Strings.EMPTY);

        Stringency = ValidationStringency.valueOf(cmd.getOptionValue(VALIDATION_STRINGENCY, DEFAULT_STRINGENCY.toString()));
    }

    public String primaryReference()
    {
        return ReferenceIds.get(0);
    }

    public List<String> allSamples()
    {
        List<String> samples = Lists.newArrayList();
        samples.addAll(ReferenceIds);
        samples.add(TumorId);
        return samples;
    }

    public int typicalReadDepth()
    {
        return DEFAULT_TYPICAL_READ_DEPTH;
    }
    public int minPartition()
    {
        return DEFAULT_MIN_PARTITION;
    }

    public boolean isValid()
    {
        if(ReferenceIds.size() != ReferenceBamPath.size())
        {
            AMB_LOGGER.error("Each reference sample must have matching bam");
            return false;
        }

        checkCreateOutputDir(OutputDir);

        if(!new File(BafLociPath).exists())
        {
            AMB_LOGGER.error("Unable to locate vcf file {}", BafLociPath);
            return false;
        }

        if(!new File(TumorBamPath).exists())
        {
            AMB_LOGGER.error("Unable to locate tumor bam file {}", TumorBamPath);
            return false;
        }

        if(!TumorOnly)
        {
            for(String referenceBam : ReferenceBamPath)
            {
                if(!new File(referenceBam).exists())
                {
                    AMB_LOGGER.error("Unable to locate reference bam {}", referenceBam);
                    return false;
                }
            }
        }

        return true;
    }

    public static Options createOptions()
    {
        final Options options = new org.apache.commons.cli.Options();
        options.addOption(TUMOR_ONLY, false, "Tumor only mode");
        options.addOption(THREADS, true, "Number of threads [" + DEFAULT_THREADS + "]");
        options.addOption(REFERENCE, true, "Name of reference sample");
        options.addOption(REFERENCE_BAM, true, "Path to reference bam file");
        options.addOption(TUMOR, true, "Name of tumor sample");
        options.addOption(TUMOR_BAM, true, "Path to tumor bam file");
        options.addOption(OUTPUT_DIR, true, "Output directory");
        options.addOption(BAF_LOCI, true, "Path to BAF loci vcf file");
        options.addOption(REF_GENOME, true, "Path to the ref genome fasta file");
        options.addOption(MIN_BASE_QUALITY, true, "Minimum quality for a base to be considered [" + DEFAULT_MIN_BASE_QUALITY + "]");
        options.addOption(MIN_MAPPING_QUALITY,
                true,
                "Minimum mapping quality for an alignment to be used [" + DEFAULT_MIN_MAPPING_QUALITY + "]");
        options.addOption(MIN_HET_AF_PERCENTAGE, true, "Min heterozygous AF% [" + DEFAULT_MIN_HET_AF_PERCENTAGE + "]");
        options.addOption(MAX_HET_AF_PERCENTAGE, true, "Max heterozygous AF% [" + DEFAULT_MAX_HET_AF_PERCENTAGE + "]");
        options.addOption(MIN_DEPTH_PERCENTAGE, true, "Max percentage of median depth [" + DEFAULT_MIN_DEPTH_PERCENTAGE + "]");
        options.addOption(MAX_DEPTH_PERCENTAGE, true, "Min percentage of median depth [" + DEFAULT_MAX_DEPTH_PERCENTAGE + "]");

        options.addOption(TUMOR_ONLY_MIN_VAF, true, "Min support in ref and alt in tumor only mode [" + DEFAULT_TUMOR_ONLY_MIN_VAF + "]");
        options.addOption(TUMOR_ONLY_MIN_SUPPORT,
                true,
                "Min VAF in ref and alt in tumor only mode [" + DEFAULT_TUMOR_ONLY_MIN_SUPPORT + "]");
        options.addOption(VALIDATION_STRINGENCY, true, "SAM validation strategy: STRICT, SILENT, LENIENT [STRICT]");

        return options;
    }

    @NotNull
    static String parameter(@NotNull final CommandLine cmd, @NotNull final String parameter, @NotNull final StringJoiner missing)
    {
        final String value = cmd.getOptionValue(parameter);
        if(value == null)
        {
            missing.add(parameter);
            return "";
        }
        return value;
    }
}
