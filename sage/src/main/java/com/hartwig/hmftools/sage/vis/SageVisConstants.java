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
}
