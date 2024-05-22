package com.hartwig.hmftools.esvee.assembly.output;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.VCF_ZIP_EXTENSION;
import static com.hartwig.hmftools.esvee.common.FileCommon.RAW_VCF_SUFFIX;

import java.util.Arrays;
import java.util.List;

import com.google.common.collect.Lists;

public enum WriteType
{
    ASSEMBLY_BAM("assembly.bam"),
    JUNC_ASSEMBLY("assemblies.tsv"),
    ASSEMBLY_READ("assembly_read.tsv"),
    BREAKEND("breakend.tsv"),
    VCF(RAW_VCF_SUFFIX),
    ALIGNMENT("alignment.tsv"),
    ALIGNMENT_DATA("align_detailed.tsv"),
    DECOY_MATCHES("decoy_matches.tsv"),
    PHASE_GROUP_BUILDING("phase_group_building.tsv");

    private final String mFileId;

    WriteType(final String fileId)
    {
        mFileId = fileId;
    }

    public String fileId() { return mFileId; }

    private static final String ALL = "ALL";
    private static final String ASSEMBLIES_STR = "ASSEMBLIES"; // for backwards compatibility

    public static List<WriteType> fromConfig(final String configStr)
    {
        List<WriteType> writeTypes = Lists.newArrayList();

        if(configStr != null)
        {
            if(configStr.equals(WriteType.ALL))
            {
                Arrays.stream(WriteType.values()).filter(x -> x != ASSEMBLY_READ).forEach(x -> writeTypes.add(x));
            }
            else
            {
                String[] writeTypeValues = configStr.split(ITEM_DELIM, -1);

                for(String writeType : writeTypeValues)
                {
                    if(writeType.equals(ASSEMBLIES_STR))
                        writeTypes.add(JUNC_ASSEMBLY);
                    else
                        writeTypes.add(WriteType.valueOf(writeType));
                }
            }
        }
        else
        {
            writeTypes.add(VCF);
            // writeTypes.add(ASSEMBLY_BAM);
        }

        return writeTypes;
    }
}
