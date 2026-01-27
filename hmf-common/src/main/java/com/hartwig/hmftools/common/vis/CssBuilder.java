package com.hartwig.hmftools.common.vis;

import static java.lang.String.format;

import java.awt.Color;
import java.util.Map;
import java.util.StringJoiner;

import org.pcollections.PMap;
import org.pcollections.TreePMap;

public final class CssBuilder
{
    public static final CssBuilder EMPTY = new CssBuilder();

    private static final String CSS_DELIM = ";";
    private static final String PROPERTY_VALUE_DELIM = ":";

    private final PMap<String, String> mProperties;
    private final PMap<Integer, SizeColorPair> mBoxShadowValues;

    private CssBuilder()
    {
        mProperties = TreePMap.empty();
        mBoxShadowValues = TreePMap.empty();
    }

    private CssBuilder(final PMap<String, String> properties, final PMap<Integer, SizeColorPair> boxShadowValues)
    {
        mProperties = properties;
        mBoxShadowValues = boxShadowValues;
    }

    private static String colorToHexString(final Color color)
    {
        return format("#%02x%02x%02x", color.getRed(), color.getGreen(), color.getBlue());
    }

    public CssBuilder merge(final CssBuilder other)
    {
        return new CssBuilder(mProperties.plusAll(other.mProperties), mBoxShadowValues.plusAll(other.mBoxShadowValues));
    }

    public CssBuilder overflow(final String value)
    {
        return new CssBuilder(mProperties.plus("overflow", value), mBoxShadowValues);
    }

    public CssBuilder textAlign(final String value)
    {
        return new CssBuilder(mProperties.plus("text-align", value), mBoxShadowValues);
    }

    public CssBuilder tableLayout(final String value)
    {
        return new CssBuilder(mProperties.plus("table-layout", value), mBoxShadowValues);
    }

    public CssBuilder borderCollapse(final String value)
    {
        return new CssBuilder(mProperties.plus("border-collapse", value), mBoxShadowValues);
    }

    public CssBuilder floatStyle(final String value)
    {
        return new CssBuilder(mProperties.plus("float", value), mBoxShadowValues);
    }

    public CssBuilder marginRight(final CssSize size)
    {
        return new CssBuilder(mProperties.plus("margin-right", size.toString()), mBoxShadowValues);
    }

    public CssBuilder writingMode(final String value)
    {
        return new CssBuilder(mProperties.plus("writing-mode", value), mBoxShadowValues);
    }

    public CssBuilder transform(final String value)
    {
        return new CssBuilder(mProperties.plus("transform", value), mBoxShadowValues);
    }

    public CssBuilder fontFamily(final String value)
    {
        return new CssBuilder(mProperties.plus("font-family", value), mBoxShadowValues);
    }

    public CssBuilder display(final String value)
    {
        return new CssBuilder(mProperties.plus("display", value), mBoxShadowValues);
    }

    public CssBuilder fontWeight(final String value)
    {
        return new CssBuilder(mProperties.plus("font-weight", value), mBoxShadowValues);
    }

    public CssBuilder width(final CssSize size)
    {
        return new CssBuilder(mProperties.plus("width", size.toString()), mBoxShadowValues);
    }

    public CssBuilder height(final CssSize size)
    {
        return new CssBuilder(mProperties.plus("height", size.toString()), mBoxShadowValues);
    }

    public CssBuilder padding(final CssSize size)
    {
        return new CssBuilder(mProperties.plus("padding", size.toString()), mBoxShadowValues);
    }

    public CssBuilder paddingRight(final CssSize size)
    {
        return new CssBuilder(mProperties.plus("padding-right", size.toString()), mBoxShadowValues);
    }

    public CssBuilder paddingLeft(final CssSize size)
    {
        return new CssBuilder(mProperties.plus("padding-left", size.toString()), mBoxShadowValues);
    }

    public CssBuilder margin(final CssSize size)
    {
        return new CssBuilder(mProperties.plus("margin", size.toString()), mBoxShadowValues);
    }

    public CssBuilder borderSpacing(final CssSize size)
    {
        return new CssBuilder(mProperties.plus("border-spacing", size.toString()), mBoxShadowValues);
    }

    public CssBuilder backgroundColor(final Color color)
    {
        return new CssBuilder(mProperties.plus("background-color", colorToHexString(color)), mBoxShadowValues);
    }

    public CssBuilder color(final Color color)
    {
        return new CssBuilder(mProperties.plus("color", colorToHexString(color)), mBoxShadowValues);
    }

    public CssBuilder borderBottom(final CssSize size, final String style, final Color color)
    {
        return new CssBuilder(mProperties.plus("border-bottom", format("%s %s %s", size.toString(), style, colorToHexString(color))), mBoxShadowValues);
    }

    public CssBuilder borderTop(final CssSize size, final String style, final Color color)
    {
        return new CssBuilder(mProperties.plus("border-top", format("%s %s %s", size.toString(), style, colorToHexString(color))), mBoxShadowValues);
    }

    public CssBuilder borderLeft(final CssSize size, final String style, final Color color)
    {
        return new CssBuilder(mProperties.plus("border-left", format("%s %s %s", size.toString(), style, colorToHexString(color))), mBoxShadowValues);
    }

    public CssBuilder borderRight(final CssSize size, final String style, final Color color)
    {
        return new CssBuilder(mProperties.plus("border-right", format("%s %s %s", size.toString(), style, colorToHexString(color))), mBoxShadowValues);
    }

    public CssBuilder border(final CssSize size, final String style, final Color color)
    {
        return new CssBuilder(mProperties.plus("border", format("%s %s %s", size.toString(), style, colorToHexString(color))), mBoxShadowValues);
    }

    public CssBuilder noBorder()
    {
        return new CssBuilder(mProperties.plus("border", "none"), mBoxShadowValues);
    }

    public CssBuilder verticalAlign(final String value)
    {
        return new CssBuilder(mProperties.plus("vertical-align", value), mBoxShadowValues);
    }

    public CssBuilder fontSizePt(int ptSize)
    {
        return new CssBuilder(mProperties.plus("font-size", ptSize + "pt"), mBoxShadowValues);
    }

    public CssBuilder fontSizePct(double pctSize)
    {
        return new CssBuilder(mProperties.plus("font-size", format("%.2f%%", 100 * pctSize)), mBoxShadowValues);
    }

    public CssBuilder leftBoxShadow(final CssSize size, final Color color)
    {
        return new CssBuilder(mProperties, mBoxShadowValues.plus(BorderDirections.LEFT.ordinal(), new SizeColorPair(size, color)));
    }

    public CssBuilder rightBoxShadow(final CssSize size, final Color color)
    {
        return new CssBuilder(mProperties, mBoxShadowValues.plus(BorderDirections.RIGHT.ordinal(), new SizeColorPair(size, color)));
    }

    private String boxShadowValuesToString()
    {
        StringJoiner joiner = new StringJoiner(", ");
        for(BorderDirections direction : BorderDirections.values())
        {
            SizeColorPair boxShadowValue = mBoxShadowValues.get(direction.ordinal());
            if(boxShadowValue == null)
            {
                continue;
            }

            CssSize hOffset = boxShadowValue.size.scale(direction.HorizontalScale);
            CssSize vOffset = boxShadowValue.size.scale(direction.VerticalScale);
            joiner.add(format("%s %s 0 0 %s", hOffset.toString(), vOffset.toString(), colorToHexString(boxShadowValue.color)));
        }

        return joiner.toString();
    }

    @Override
    public String toString()
    {
        if(mProperties.isEmpty() && mBoxShadowValues.isEmpty())
        {
            return "";
        }

        StringJoiner styleJoiner = new StringJoiner(CSS_DELIM, "", CSS_DELIM);
        for(Map.Entry<String, String> propertyEntry : mProperties.entrySet())
        {
            String property = propertyEntry.getKey();
            String value = propertyEntry.getValue();
            styleJoiner.add(property + PROPERTY_VALUE_DELIM + value);
        }

        if(!mBoxShadowValues.isEmpty())
        {
            String boxShadowValue = boxShadowValuesToString();
            styleJoiner.add("box-shadow" + PROPERTY_VALUE_DELIM + boxShadowValue);
        }

        return styleJoiner.toString();
    }

    private enum BorderDirections
    {
        TOP(0, -1),
        RIGHT(1, 0),
        BOTTOM(0, 1),
        LEFT(-1, 0);

        public int HorizontalScale;
        public int VerticalScale;

        BorderDirections(int horizontalScale, int verticalScale)
        {
            HorizontalScale = horizontalScale;
            VerticalScale = verticalScale;
        }
    }

    private static class SizeColorPair
    {
        public final CssSize size;
        public final Color color;

        public SizeColorPair(final CssSize size, final Color color)
        {
            this.size = size;
            this.color = color;
        }
    }
}
