package com.hartwig.hmftools.esvee;

public enum WriteType
{
    ASSEMBLY_BAM("assembly.bam"),
    ASSEMBLIES("assemblies.tsv"),
    ASSEMBLY_READS("assembly_reads.tsv"),
    READS("reads.tsv"),
    VCF("esvee.vcf.gz"),
    BREAKENDS("breakends.tsv"); // not currently defined or written

    private final String mFileId;

    WriteType(final String fileId)
    {
        mFileId = fileId;
    }

    public String fileId() { return mFileId; }

    public static final String ALL = "ALL";
}
