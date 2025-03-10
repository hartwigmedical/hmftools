package com.hartwig.hmftools.linx;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;

public enum WriteType
{
    SV_DATA,
    CLUSTER,
    LINK,
    ECDNA,
    CLUSTERING,
    VIS_DATA;

    public static final String WRITE_ALL = "ALL";
    public static final String WRITE_STANDARD = "STANDARD";

    public static final List<WriteType> STANDARD_TYPES = Lists.newArrayList(SV_DATA, CLUSTER, LINK);

    public static List<WriteType> parseConfig(final String configStr)
    {
        if(configStr.equals(WRITE_ALL))
            return Arrays.stream(WriteType.values()).collect(Collectors.toList());

        if(configStr.equals(WRITE_STANDARD))
            return STANDARD_TYPES;

        String[] configItems = configStr.split(ITEM_DELIM, -1);
        List<WriteType> writeTypes = Lists.newArrayList();

        for(String configItem : configItems)
        {
            writeTypes.add(WriteType.valueOf(configItem));
        }

        return writeTypes;
    }
}
