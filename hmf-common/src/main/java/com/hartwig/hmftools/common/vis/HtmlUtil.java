package com.hartwig.hmftools.common.vis;

import static j2html.TagCreator.table;

import java.util.List;

import j2html.tags.DomContent;

public final class HtmlUtil
{
    // sizes
    private static final int BASE_FONT_SIZE = 10;

    // styles
    public static final CssBuilder BASE_FONT_STYLE = CssBuilder.EMPTY.fontSizePt(BASE_FONT_SIZE).fontFamily("sans-serif");

    private HtmlUtil() {}


    public static DomContent styledTable(final List<DomContent> elems, final CssBuilder style)
    {
        return table().with(elems).withStyle(BASE_FONT_STYLE.merge(style).toString());
    }
}
