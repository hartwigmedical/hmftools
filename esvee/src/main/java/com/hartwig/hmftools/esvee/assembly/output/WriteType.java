package com.hartwig.hmftools.esvee.assembly.output;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.esvee.common.FileCommon.RAW_VCF_SUFFIX;

import java.util.Arrays;
import java.util.List;

import com.google.common.collect.Lists;

public enum WriteType
{
    ASSEMBLY_BAM("assembly.bam"),
    JUNC_ASSEMBLY("assembly.tsv"),
    ASSEMBLY_READ("assembly_read.tsv"),
    BREAKEND("breakend.tsv"),
    VCF(RAW_VCF_SUFFIX),
    PHASED_ASSEMBLY("phased_assembly.tsv"),
    ALIGNMENT("alignment.tsv"),
    DECOY_MATCHES("decoy_match.tsv"),
    PHASE_GROUP_BUILDING("phase_group_building.tsv");

    private final String mFileId;

    WriteType(final String fileId)
    {
        mFileId = fileId;
    }

    public String fileId() { return mFileId; }

    public static boolean requiresAlignment(final List<WriteType> writeTypes)
    {
        for(WriteType writeType : writeTypes)
        {
            switch(writeType)
            {
                case BREAKEND:
                case ALIGNMENT:
                case VCF:
                    return true;

                default:
                    break;
            }
        }

        return false;
    }

    private static final String ALL = "ALL";

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
                    // backwards compatibility
                    if(writeType.equals("ALIGNMENTS"))
                        writeTypes.add(ALIGNMENT);
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
