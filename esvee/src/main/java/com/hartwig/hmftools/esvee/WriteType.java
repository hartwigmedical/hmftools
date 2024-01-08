package com.hartwig.hmftools.esvee;

public enum WriteType
{
    ASSEMBLY_BAM("assembly.bam"),
    ASSEMBLIES("assemblies.tsv"),
    READS("reads.tsv"),
    BREAKENDS("breakends.tsv");

    private final String mFileId;

    WriteType(final String fileId)
    {
        mFileId = fileId;
    }

    public String fileId() { return mFileId; }

    public static final String ALL = "ALL";
}
