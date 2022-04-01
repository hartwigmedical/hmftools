package com.hartwig.hmftools.amber;

import static com.hartwig.hmftools.amber.AmberConstants.*;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkCreateOutputDir;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import com.beust.jcommander.IStringConverter;
import com.google.common.collect.Lists;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import com.beust.jcommander.Parameter;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import htsjdk.samtools.ValidationStringency;

public class AmberConfig
{
    public static final int DEFAULT_THREADS = 1;

    @Parameter(names = "-tumor", description = "Name of tumor sample")
    public String TumorId;

    @Parameter(names = "-tumor_bam", description = "Path to indexed tumor bam/cram file")
    public String TumorBamPath;

    @Parameter(names = "-reference", description = "Name of reference sample")
    public List<String> ReferenceIds = new ArrayList<>();

    @Parameter(names = "-reference_bam", description = "Path to reference bam/cram file")
    public List<String> ReferenceBamPath = new ArrayList<>();

    @Parameter(names = "-loci", required = true, description = "Path to BAF loci vcf file")
    public String BafLociPath;

    @Parameter(names = "-" + REF_GENOME,
               description = "Path to the reference genome fasta file. Required only when using CRAM files.")
    public String RefGenomePath;

    @Parameter(names = "-output_dir",
               required = true,
               description = "Path to the output directory. "
                       + "This directory will be created if it does not already exist.")
    public String OutputDir;

    @Parameter(names = "-tumor_only_min_support",
               description = "Min support in ref and alt in tumor only mode")
    public int TumorOnlyMinSupport = DEFAULT_TUMOR_ONLY_MIN_SUPPORT;

    @Parameter(names = "-tumor_only_min_vaf",
               description = "Min VAF in ref and alt in tumor only mode")
    public double TumorOnlyMinVaf = DEFAULT_TUMOR_ONLY_MIN_VAF;

    @Parameter(names = "-tumor_only_min_depth",
               description = "Min depth in tumor only mode")
    public int TumorOnlyMinDepth = DEFAULT_TUMOR_ONLY_MIN_DEPTH;

    @Parameter(names = "-min_base_quality",
               description = "Minimum quality for a base to be considered")
    public int MinBaseQuality = DEFAULT_MIN_BASE_QUALITY;

    @Parameter(names = "-min_mapping_quality",
               description = "Minimum mapping quality for an alignment to be used")
    public int MinMappingQuality = DEFAULT_MIN_MAPPING_QUALITY;

    @Parameter(names = "-min_depth_percent",
               description = "Min percentage of median depth")
    public double MinDepthPercent = DEFAULT_MIN_DEPTH_PERCENTAGE;

    @Parameter(names = "-max_depth_percent",
               description = "Max percentage of median depth")
    public double MaxDepthPercent = DEFAULT_MAX_DEPTH_PERCENTAGE;

    @Parameter(names = "-min_het_af_percent",
               description = "Min heterozygous AF%")
    public double MinHetAfPercent = DEFAULT_MIN_HET_AF_PERCENTAGE;

    @Parameter(names = "-max_het_af_percent",
               description = "Max heterozygous AF%")
    public double MaxHetAfPercent = DEFAULT_MAX_HET_AF_PERCENTAGE;

    @Parameter(names = "-validation_stringency",
               description = "SAM validation strategy")
    public ValidationStringency Stringency = ValidationStringency.DEFAULT_STRINGENCY;

    @Parameter(names = "-threads",
               description = "Number of threads")
    public int ThreadCount = DEFAULT_THREADS;

    @Parameter(names = "-" + RefGenomeVersion.REF_GENOME_VERSION,
               required = true,
               description = RefGenomeVersion.REF_GENOME_VERSION_CFG_DESC,
               converter = RefGenomeVersionConverter.class)
    public RefGenomeVersion refGenomeVersion;

    public static final Logger AMB_LOGGER = LogManager.getLogger(AmberConfig.class);

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

    public boolean isTumorOnly() { return ReferenceBamPath.isEmpty() && TumorBamPath != null; }

    public boolean isGermlineOnly()
    {
        return !ReferenceBamPath.isEmpty() && TumorBamPath == null;
    }

    // use the tumor id if it is not null, otherwise primary reference Id
    public String getSampleId()
    {
        return TumorId != null ? TumorId : primaryReference();
    }

    public boolean isValid()
    {
        if(ReferenceIds.size() != ReferenceBamPath.size())
        {
            AMB_LOGGER.error("Each reference sample must have matching bam");
            return false;
        }

        if ((TumorId == null) != (TumorBamPath == null))
        {
            AMB_LOGGER.error("Unmatched: TumorId: {} and TumorBamPath: {}", TumorId, TumorBamPath);
            return false;
        }

        checkCreateOutputDir(OutputDir);

        if(!new File(BafLociPath).exists())
        {
            AMB_LOGGER.error("Unable to locate vcf file {}", BafLociPath);
            return false;
        }

        if(TumorBamPath != null && !new File(TumorBamPath).exists())
        {
            AMB_LOGGER.error("Unable to locate tumor bam file {}", TumorBamPath);
            return false;
        }

        if(!isTumorOnly())
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

    // we need to define a converter for ref genome version
    static class RefGenomeVersionConverter implements IStringConverter<RefGenomeVersion>
    {
        @Override
        public RefGenomeVersion convert(String value)
        {
            return RefGenomeVersion.from(value);
        }
    }
}
