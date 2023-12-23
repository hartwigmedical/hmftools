package com.hartwig.hmftools.esvee.output.html;

import java.util.ArrayList;
import java.util.List;

import org.jetbrains.annotations.Nullable;

public class PieChart
{
    private static final int DEFAULT_PIE_SCALE = 150;
    private static final int DEFAULT_BORDER_BUFFER = 20;
    private static final int DEFAULT_FONT_SIZE = 16;
    private static final List<String> COLOURS = List.of("red", "green", "blue", "yellow", "purple", "orange", "lime",
            "hotpink", "aqua", "brown", "yellowgreen", "pink", "lightblue");

    private static class Slice
    {
        public final String Name;
        public final long Amount;
        public final String AmountLabel;

        private Slice(final String name, final long amount, final String amountLabel)
        {
            Name = name;
            Amount = amount;
            AmountLabel = amountLabel;
        }
    }

    private final String mTitle;
    private final double mPieScale;
    private final double mBorderBuffer;
    private final double mFontSize;
    private final double mMinSliceWidthToLabel = Math.toRadians(10);
    private final List<Slice> mSlices = new ArrayList<>();

    public PieChart(final String title)
    {
        this(title, DEFAULT_PIE_SCALE, DEFAULT_BORDER_BUFFER, DEFAULT_FONT_SIZE);
    }

    public PieChart(final String title, final double pieScale, final double borderBuffer, final double fontSize)
    {
        mTitle = title;
        mPieScale = pieScale;
        mBorderBuffer = borderBuffer;
        mFontSize = fontSize;
    }

    public void add(final String item, final long amount)
    {
        add(item, amount, String.valueOf(amount));
    }

    public void add(final String item, final long amount, final String amountLabel)
    {
        if (amount < 0)
            throw new IllegalArgumentException("A pie slice cannot represent a negative amount");
        else if (amount == 0)
            return;

        mSlices.add(new Slice(item, amount, amountLabel));
    }

    public String toSVG()
    {
        return appendAsSVG(new HTMLBuilder(), "").toString();
    }

    public HTMLBuilder appendAsSVG(final HTMLBuilder builder, final String userStyle)
    {
        final int legendWidth = calculateLegendWidth();
        final int legendHeight = (int) (mSlices.size() * mFontSize * 2);
        final double titleHeight = mFontSize * 2;

        final double xOffset = -mPieScale - mBorderBuffer;
        final double width = mPieScale * 2 + mBorderBuffer * 3 + legendWidth;
        builder.appendStartTag("<svg class=\"pieChart\" style=\"%s\" viewBox=\"%s %s %s %s\">", userStyle,
                xOffset, -mPieScale - mBorderBuffer - titleHeight,
                width,
                Math.max(mPieScale * 2 + mBorderBuffer * 2, mBorderBuffer + legendHeight) + titleHeight);
        builder.appendStartTag("<style>");
        builder.append(".pieChart text {\n");
        builder.append("\tfill: black;\n");
        builder.append("\ttext-anchor: middle;\n");
        builder.append("\tdominant-baseline: central;\n");
        builder.append("\tfont-size: %spx;\n", mFontSize);
        builder.append("}\n");
        builder.append(".pieChart text.title {\n");
        builder.append("\tfont-size: %spx;\n", mFontSize * 2);
        builder.append("}\n");
        builder.append(".pieChart .pieSlice path {\n");
        builder.append("\tstroke: black;\n");
        builder.append("\tstroke-width: 0.5px;\n");
        builder.append("}\n");
        builder.append(".pieChart .legend {\n");
        builder.append("\tfont-size: %spx;\n", mFontSize);
        builder.append("\tfont-family: monospace;\n");
        builder.append("}\n");
        builder.appendEndTag("</style>");

        builder.appendStartTag("<g>");
        final double titleX = width / 2 + xOffset;
        final double titleY = -mPieScale - mBorderBuffer - (titleHeight / 2);
        builder.append("<text class=\"title\" transform=\"translate(%s, %s)\">%s</text>\n", titleX, titleY, mTitle);
        builder.appendEndTag("<g>");

        builder.appendStartTag("<g>");

        final long total = mSlices.stream().mapToLong(slice -> slice.Amount).sum();
        final var labelBuilder = new HTMLBuilder();
        final var legendBuilder = new HTMLBuilder();
        legendBuilder.appendStartTag("<table class=\"legend\">");
        if (total > 0)
        {
            double fractionSoFar = 0;
            for(int i = 0; i < mSlices.size(); i++)
            {
                final Slice slice = mSlices.get(i);
                final boolean isLastSlice = i + 1 == mSlices.size();

                builder.appendStartTag("<g class=\"pieSlice\">");
                builder.append("<title>").append(slice.Name).append(": ").append(slice.AmountLabel).append("</title>\n");
                final double fractionOfTotal = ((double) slice.Amount) / total;
                final double startFraction = fractionSoFar;
                final double endFraction = isLastSlice ? 1 : startFraction + fractionOfTotal;
                final String colour = COLOURS.get(i % COLOURS.size());
                appendSegment(builder, labelBuilder, startFraction, endFraction, colour, slice.AmountLabel);
                builder.appendEndTag("</g>");

                final String block = "&#x2588;";
                legendBuilder.append("<tr><td style=\"color: %s\">%s</td><td>%s</td><td>%s</td></tr>\n",
                        colour, block, slice.Name, slice.AmountLabel);

                fractionSoFar = endFraction;
            }
        }
        builder.appendStartTag("<g class=\"pieLabels\">");
        builder.append(labelBuilder);
        builder.appendEndTag("</g>");

        builder.appendStartTag("<g class=\"pieLegend\">");
        legendBuilder.appendEndTag("</table>");
        builder.appendStartTag("<foreignObject x=\"%s\" y=\"%s\" width=\"%s\" height=\"%s\">",
                mPieScale + mBorderBuffer, -mPieScale,
                legendWidth, legendHeight);
        builder.append(legendBuilder.toString());
        builder.appendEndTag("</foreignObject>");
        builder.appendEndTag("</g>");

        builder.appendEndTag("</g>");

        builder.appendEndTag("</svg>");
        return builder;
    }

    private int calculateLegendWidth()
    {
        final int amountLabelWidth = (int) mSlices.stream().mapToDouble(slice -> slice.AmountLabel.length() * mFontSize / 1.5).max().orElse(0);
        final int nameLabelWidth = (int) mSlices.stream().mapToDouble(slice -> slice.Name.length() * mFontSize / 1.5).max().orElse(0);
        return nameLabelWidth + amountLabelWidth;
    }

    private void appendSegment(final HTMLBuilder pathBuilder, final HTMLBuilder labelBuilder,
            final double startFraction, final double endFraction, final String fill, @Nullable final String label)
    {
        final double startRadians = 2 * Math.PI * startFraction;
        final double endRadians = 2 * Math.PI * endFraction;

        final double startX = Math.sin(startRadians) * mPieScale;
        final double startY = -Math.cos(startRadians) * mPieScale;

        final double endX = Math.sin(endRadians) * mPieScale;
        final double endY = -Math.cos(endRadians) * mPieScale;

        final double xRotation = 0;
        final boolean isLargeArc = endRadians - startRadians > Math.PI;
        final boolean sweepFlag = true;

        //noinspection ConstantValue
        pathBuilder.append("<path d=\"M 0 0 L %s %s A %s %s %s %s %s %s %s Z\" fill=\"%s\"/>\n", startX, startY, mPieScale, mPieScale, xRotation,
                isLargeArc ? 1 : 0, sweepFlag ? 1 : 0, endX, endY, fill);

        if (label == null)
            return;

        if (endRadians - startRadians < mMinSliceWidthToLabel)
            return;

        final double midRadians = (endRadians - startRadians) / 2 + startRadians;
        final double labelX = Math.sin(midRadians) * (mPieScale * 0.7);
        final double labelY = -Math.cos(midRadians) * (mPieScale * 0.7);
        final double midDegrees = Math.toDegrees(midRadians);
        final double rotationDegrees = midDegrees > 180
                ? midDegrees - 270
                : midDegrees - 90;

        labelBuilder.append("<text transform=\"translate(%s, %s) rotate(%s)\">%s</text>\n", labelX, labelY, rotationDegrees, label);
    }
}
