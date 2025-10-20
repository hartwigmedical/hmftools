package com.hartwig.hmftools.sage.vis;

import java.util.EnumMap;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.sage.common.ReadContextMatch;

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
    }

    public static final String IMPACT_KEY = "IMPACT";
    public static final int GENE_NAME_IDX = 0;
    public static final int TRANSCRIPT_NAME_IDX = 1;
    public static final int AA_VARIANT_TYPE_IDX = 2;
    public static final int HGVS_INDEX = 6;
}
