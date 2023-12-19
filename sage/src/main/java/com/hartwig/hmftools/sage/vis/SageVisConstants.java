package com.hartwig.hmftools.sage.vis;

import java.util.EnumMap;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;

public class SageVisConstants
{
    // sizes
    public static final double READ_HEIGHT_PX = 12.0;
    public static final CssSize VARIANT_INFO_SPACING_SIZE = CssSize.em(2);
    public static final int BASE_FONT_SIZE = 10;

    // styles
    public static final CssBuilder BASE_FONT_STYLE = CssBuilder.EMPTY.fontSizePt(BASE_FONT_SIZE).fontFamily("sans-serif");

    // config
    public static final int READ_EXTEND_LENGTH = 75;
    public static final int DISPLAY_EVERY_NTH_COORD = 10;
    public static final int MAX_MAPQ_SHADING_CUTTOFF = 60;
    public static final int INSERT_SIZE_CUTTOFF = 1000;

    public static final int MAX_READ_UPPER_LIMIT = 1000;
    public static final EnumMap<ReadContextCounter.MatchType, Integer> MAX_READS_PER_TYPE;

    static
    {
        MAX_READS_PER_TYPE = Maps.newEnumMap(ReadContextCounter.MatchType.class);
        MAX_READS_PER_TYPE.put(ReadContextCounter.MatchType.NONE, 5);
        MAX_READS_PER_TYPE.put(ReadContextCounter.MatchType.FULL, 40);
        MAX_READS_PER_TYPE.put(ReadContextCounter.MatchType.PARTIAL, 10);
        MAX_READS_PER_TYPE.put(ReadContextCounter.MatchType.CORE, 10);
        MAX_READS_PER_TYPE.put(ReadContextCounter.MatchType.REALIGNED, 20);
        MAX_READS_PER_TYPE.put(ReadContextCounter.MatchType.CORE_PARTIAL, 10);
        MAX_READS_PER_TYPE.put(ReadContextCounter.MatchType.REF, 20);
        MAX_READS_PER_TYPE.put(ReadContextCounter.MatchType.ALT, 5);
    }
}
