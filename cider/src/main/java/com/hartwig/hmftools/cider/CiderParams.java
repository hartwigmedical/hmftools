package com.hartwig.hmftools.cider;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;

import com.beust.jcommander.IStringConverter;
import com.beust.jcommander.Parameter;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import htsjdk.samtools.ValidationStringency;

public class CiderParams
{
    public static final int DEFAULT_THREADS = 1;
    public static final int DEFAULT_MAX_ANCHOR_ALIGN_DISTANCE = 150;

    @Parameter(names = "-id", description = "Name of sample")
    public String SampleId;

    @Parameter(names = "-bam", description = "Path to indexed bam/cram file")
    public String BamPath;

    @Parameter(names = "-" + REF_GENOME,
               description = "Path to the reference genome fasta file. Required only when using CRAM files.")
    public String RefGenomePath;

    @Parameter(names = "-output_dir",
               required = true,
               description = "Path to the output directory. "
                       + "This directory will be created if it does not already exist.")
    public String OutputDir;

    @Parameter(names = "-validation_stringency",
               description = "SAM validation strategy")
    public ValidationStringency Stringency = ValidationStringency.DEFAULT_STRINGENCY;

    @Parameter(names = "-threads",
               description = "Number of threads")
    public int ThreadCount = DEFAULT_THREADS;

    @Parameter(names = "-anchor_align_distance",
               description = "Anchor align distance in the gene location")
    public int MaxAnchorAlignDistance = DEFAULT_MAX_ANCHOR_ALIGN_DISTANCE;

    @Parameter(names = "-" + RefGenomeVersion.REF_GENOME_VERSION,
               required = true,
               description = RefGenomeVersion.REF_GENOME_VERSION_CFG_DESC,
               converter = RefGenomeVersionConverter.class)
    public RefGenomeVersion refGenomeVersion;

    @Parameter(names = "-ensembl_data_dir",
               required = true,
               description = "Ensembl data file directory")
    public String ensemblDataDir;

    @Parameter(names = "-min_base_quality",
               description = "Minimum quality for a base to be considered")
    public int MinBaseQuality = 25;

    @Parameter(names = "-assume_read_pair_overlap",
               description = "Assume that reads of the same fragment overlap. This can be used for pairing")
    public boolean assumeReadPairOverlap = false;

    @Parameter(names = "-specific_chr",
               description = "Limit to a specific chromosome")
    public String SpecificChr;

    @Parameter(names = "-write_filtered_bam",
               description = "Write a output BAM file containing all CDR3 reads")
    public boolean writeFilteredBam;

    @Parameter(names = "-num_trim_bases",
               description = "Number of bases to trim on each side of reads")
    public int numBasesToTrim = 0;

    public boolean isValid()
    {
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
