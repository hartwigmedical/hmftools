package com.hartwig.hmftools.esvee.output;

import java.util.List;

public enum WriteType
{
    ASSEMBLY_BAM("assembly.bam"),
    ASSEMBLIES("assemblies.tsv"),
    READS("assembly_reads.tsv"),
    VCF("esvee.vcf.gz"),
    DECOY_MATCHES("decoy_matches.tsv"),
    PHASE_GROUP_BUILDING("phase_group_building.tsv");

    private final String mFileId;

    WriteType(final String fileId)
    {
        mFileId = fileId;
    }

    public String fileId() { return mFileId; }

    public static final String ALL = "ALL";

    public static boolean hasAssemblyOutput(final List<WriteType> writeTypes)
    {
        return writeTypes.contains(ASSEMBLIES) || writeTypes.contains(READS) || writeTypes.contains(ASSEMBLY_BAM);
    }
}
