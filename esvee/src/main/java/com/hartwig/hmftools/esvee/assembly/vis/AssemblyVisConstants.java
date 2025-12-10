package com.hartwig.hmftools.esvee.assembly.vis;

import com.hartwig.hmftools.common.vis.CssBuilder;

public final class AssemblyVisConstants
{
    private AssemblyVisConstants() {}

    // sizes
    public static final double READ_HEIGHT_PX = 12.0;
    private static final int BASE_FONT_SIZE = 10;

    // styles
    public static final CssBuilder BASE_FONT_STYLE = CssBuilder.EMPTY.fontSizePt(BASE_FONT_SIZE).fontFamily("sans-serif");
}
