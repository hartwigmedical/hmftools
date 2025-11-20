package com.hartwig.hmftools.esvee.common;

import static com.hartwig.hmftools.common.utils.config.ConfigItem.enumValueSelectionAsStr;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.esvee.common.FileCommon.RAW_VCF_SUFFIX;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.PREP_DISC_STATS_FILE_ID;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.PREP_FRAG_LENGTH_FILE_ID;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.PREP_JUNCTION_FILE_ID;

import java.util.Arrays;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public enum WriteType
{
    // prep routine
    PREP_JUNCTION(PREP_JUNCTION_FILE_ID, true),
    PREP_READ("read.tsv", true),
    PREP_BAM("bam", true),
    UNSORTED_BAM("unsorted.bam", true),
    CACHE_BAM("cache", true),
    FRAGMENT_LENGTH_DIST(PREP_FRAG_LENGTH_FILE_ID, true),
    DISCORDANT_STATS(PREP_DISC_STATS_FILE_ID, true),

    // assembly routine
    ASSEMBLY_BAM("assembly.bam", false),
    JUNC_ASSEMBLY("assembly.tsv", false),
    ASSEMBLY_READ("assembly_read.tsv", false),
    BREAKEND("breakend.tsv", false),
    VCF(RAW_VCF_SUFFIX, false),
    PHASED_ASSEMBLY("phased_assembly.tsv", false),
    ALIGNMENT("alignment.tsv", false),
    DECOY_MATCHES("decoy_match.tsv", false),
    PHASE_GROUP_BUILDING("phase_group_building.tsv", false);

    private final String mFileId;
    private final boolean mPrepType; // otherwise assembly

    WriteType(final String fileId, final boolean isPrepType)
    {
        mFileId = fileId;
        mPrepType = isPrepType;
    }

    public String fileId() { return mFileId; }
    public boolean isPrepType() { return mPrepType; }
    public boolean isAssemblyType() { return !mPrepType; }

    private static final String PREP_STANDARD_CFG = "PREP_STANDARD";

    private static final String ASSEMBLY_ALL_CFG = "ASSEMBLY_ALL";
    private static final String ASSEMBLY_STANDARD_CFG = "ASSEMBLY_STANDARD";

    private static final Set<WriteType> PREP_STANDARD = Set.of(PREP_JUNCTION, PREP_BAM, FRAGMENT_LENGTH_DIST, DISCORDANT_STATS);
    private static final Set<WriteType> ASSEMBLY_STANDARD = Set.of(VCF, BREAKEND, JUNC_ASSEMBLY);

    public static Set<WriteType> parsePrepTypes(final String configStr)
    {
        Set<WriteType> writeTypes = Sets.newHashSet();

        if(configStr == null || configStr.isEmpty())
            return PREP_STANDARD;

        if(configStr.contains(PREP_STANDARD_CFG))
            writeTypes.addAll(PREP_STANDARD);

        parseTypes(writeTypes, configStr);

        return writeTypes.stream().filter(x -> x.isPrepType()).collect(Collectors.toSet());
    }

    public static Set<WriteType> parseAssemblyTypes(final String configStr)
    {
        Set<WriteType> writeTypes = Sets.newHashSet();

        if(configStr == null || configStr.isEmpty())
            return ASSEMBLY_STANDARD;

        if(configStr.contains(ASSEMBLY_STANDARD_CFG))
            writeTypes.addAll(ASSEMBLY_STANDARD);

        if(configStr.contains(ASSEMBLY_ALL_CFG))
        {
            Arrays.stream(WriteType.values()).forEach(x -> writeTypes.add(x));
        }

        parseTypes(writeTypes, configStr);

        return writeTypes.stream().filter(x -> x.isAssemblyType()).collect(Collectors.toSet());
    }

    private static void parseTypes(final Set<WriteType> writeTypes, final String configStr)
    {
        String[] writeTypesArray = configStr.split(ITEM_DELIM, -1);

        for(String writeTypeStr : writeTypesArray)
        {
            if(writeTypeStr.equals(PREP_STANDARD_CFG) || writeTypeStr.equals(ASSEMBLY_STANDARD_CFG) || writeTypeStr.equals(ASSEMBLY_ALL_CFG))
                continue;

            WriteType writeType = WriteType.valueOf(writeTypeStr);

            if(!writeTypes.contains(writeType))
                writeTypes.add(writeType);
        }
    }

    public static final String WRITE_TYPES = "write_types";

    public static void registerWriteTypes(final ConfigBuilder configBuilder)
    {
        if(!configBuilder.isRegistered(WRITE_TYPES))
            configBuilder.addConfigItem(WRITE_TYPES, false, enumValueSelectionAsStr(WriteType.values(), "Write types"));
    }
}
