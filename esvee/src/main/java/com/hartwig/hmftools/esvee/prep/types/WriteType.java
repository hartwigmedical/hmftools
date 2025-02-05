package com.hartwig.hmftools.esvee.prep.types;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;

import java.util.List;

import com.google.common.collect.Lists;

public enum WriteType
{
    PREP_JUNCTION,
    PREP_READ,
    PREP_BAM,
    UNSORTED_BAM,
    CACHE_BAM,
    FRAGMENT_LENGTH_DIST,
    DISCORDANT_STATS;

    private static final String STANDARD_TYPES = "PREP_STANDARD";

    public static List<WriteType> parseConfigStr(final String configStr)
    {
        if(configStr == null || configStr.isEmpty() || configStr.contains(STANDARD_TYPES))
            return defaultTypes();

        String[] writeTypesArray = configStr.split(ITEM_DELIM, -1);

        List<WriteType> writeTypes = Lists.newArrayList();

        for(String writeTypeStr : writeTypesArray)
        {
            try
            {
                writeTypes.add(WriteType.valueOf(writeTypeStr));
            }
            catch(Exception e) {} // invalid or may be an assembly write type
        }

        return writeTypes;
    }

    public static List<WriteType> defaultTypes() { return List.of(PREP_JUNCTION, PREP_BAM, FRAGMENT_LENGTH_DIST, DISCORDANT_STATS); }
}
