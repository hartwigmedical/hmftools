package com.hartwig.hmftools.sage.vis;

import java.util.EnumMap;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.sage.common.ReadContextMatch;

public final class SageVisConstants
{
    private SageVisConstants() {}

    // sizes
    public static final double READ_HEIGHT_PX = 12.0;

    // config
    public static final int READ_EXTEND_LENGTH = 75;
    public static final int DISPLAY_EVERY_NTH_COORD = 10;
    public static final int MAX_MAPQ_SHADING_CUTTOFF = 60;
    public static final int INSERT_SIZE_CUTTOFF = 1000;

    public static final int MAX_READ_UPPER_LIMIT = 1000;
    public static final EnumMap<ReadContextMatch, Integer> MAX_READS_PER_TYPE;

    static
    {
        MAX_READS_PER_TYPE = Maps.newEnumMap(ReadContextMatch.class);
        MAX_READS_PER_TYPE.put(ReadContextMatch.NONE, 5);
        MAX_READS_PER_TYPE.put(ReadContextMatch.FULL, 40);
        MAX_READS_PER_TYPE.put(ReadContextMatch.PARTIAL_CORE, 10);
        MAX_READS_PER_TYPE.put(ReadContextMatch.CORE, 10);
        MAX_READS_PER_TYPE.put(ReadContextMatch.REALIGNED, 20);
        MAX_READS_PER_TYPE.put(ReadContextMatch.REF, 20);
        MAX_READS_PER_TYPE.put(ReadContextMatch.SIMPLE_ALT, 5);
    }

    public static final int REF_BUFFER_SIZE = 100;
}
