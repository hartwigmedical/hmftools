package com.hartwig.hmftools.sage.vis;

import static java.lang.String.format;
import static java.util.Map.entry;

import static com.hartwig.hmftools.sage.vis.ColorUtil.PURPLE;
import static com.hartwig.hmftools.sage.vis.ColorUtil.lighten;
import static com.hartwig.hmftools.sage.vis.SvgUtil.drawStringFromCenter;
import static com.hartwig.hmftools.sage.vis.SvgUtil.getStringBounds;

import java.awt.Color;
import java.awt.Font;
import java.awt.geom.AffineTransform;
import java.awt.geom.Path2D;
import java.awt.geom.Rectangle2D;
import java.util.List;
import java.util.Locale;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.BaseRegion;

import org.jetbrains.annotations.Nullable;
import org.jfree.svg.SVGGraphics2D;

public class SvgRender
{
    public static final Color FORWARD_STRAND_COLOR = new Color(81, 144, 207);
    public static final Color REVERSE_STRAND_COLOR = new Color(242, 167, 121);
    public static final Color OVERLAPPING_FRAGMENT_BORDER_COLOR = lighten(PURPLE, 0.5f);
    private static final Color INSERT_COLOR = new Color(153, 50, 204);
    private static final Map<Character, Color> BASE_BG_COLOUR = Map.ofEntries(
            entry('G', new Color(255, 232, 145)),
            entry('C', new Color(212, 221, 240)),
            entry('A', new Color(191, 236, 199)),
            entry('T', new Color(255, 192, 199))
    );
    private static final Map<Character, Color> BASE_FG_COLOUR = Map.ofEntries(
            entry('G', new Color(151, 83, 9)),
            entry('C', new Color(11, 105, 186)),
            entry('A', new Color(15, 95, 16)),
            entry('T', new Color(145, 0, 7))
    );

    private static final int BOX_PADDING = 1;
    private static final double BASE_BOX_SIZE = 30.0;
    private static final double DEL_CONNECTOR_BOX_PROPORTION = 0.2;
    private static final double BOUNDARY_BOX_PROPORTION = 0.2;

    // font sizes are tied to BASE_BOX_SIZE
    private static final int INDEL_FONT_SIZE = 18;
    private static final int BASE_FONT_SIZE = 26;
    private static final int COORD_FONT_SIZE = 30;
    private static final Font INDEL_FONT = new Font("SansSerif", Font.PLAIN, INDEL_FONT_SIZE);
    private static final Font BASE_FONT = new Font("SansSerif", Font.PLAIN, BASE_FONT_SIZE);
    private static final Font COORD_FONT = new Font("SansSerif", Font.PLAIN, COORD_FONT_SIZE);

    private static final int MAX_BASEQ_SHADING_CUTTOFF = 37;

    private static void drawTopBoxBorder(SVGGraphics2D svgCanvas, int boxIdx, double boxPropWidth)
    {
        // canvas is scaled to be in units of BASE_BOX_SIZE
        Rectangle2D rect = new Rectangle2D.Double(boxIdx, 0.0, 1.0, boxPropWidth);
        svgCanvas.fill(rect);
    }

    private static void drawRightBoxBorder(SVGGraphics2D svgCanvas, int boxIdx, double boxPropWidth)
    {
        // canvas is scaled to be in units of BASE_BOX_SIZE
        Rectangle2D rect = new Rectangle2D.Double(boxIdx + 1 - boxPropWidth, 0.0, boxPropWidth, 1.0);
        svgCanvas.fill(rect);
    }

    private static void drawBottomBoxBorder(SVGGraphics2D svgCanvas, int boxIdx, double boxPropWidth)
    {
        // canvas is scaled to be in units of BASE_BOX_SIZE
        Rectangle2D rect = new Rectangle2D.Double(boxIdx, 1.0 - boxPropWidth, 1.0, boxPropWidth);
        svgCanvas.fill(rect);
    }

    private static void drawLeftBoxBorder(SVGGraphics2D svgCanvas, int boxIdx, double boxPropWidth)
    {
        // canvas is scaled to be in units of BASE_BOX_SIZE
        Rectangle2D rect = new Rectangle2D.Double(boxIdx, 0.0, boxPropWidth, 1.0);
        svgCanvas.fill(rect);
    }

    private static void drawRightInsertIndicator(final SVGGraphics2D svgCanvas, int boxIdx)
    {
        // canvas is scaled to be in units of BASE_BOX_SIZE
        drawRightBoxBorder(svgCanvas, boxIdx, BOUNDARY_BOX_PROPORTION);

        Rectangle2D rightTopI = new Rectangle2D.Double(boxIdx + 1 - 1.0 / 3, 0.0, 1.0 / 3, BOUNDARY_BOX_PROPORTION);
        svgCanvas.fill(rightTopI);

        Rectangle2D rightBottomI =
                new Rectangle2D.Double(boxIdx + 1 - 1.0 / 3, 1.0 - BOUNDARY_BOX_PROPORTION, 1.0 / 3, BOUNDARY_BOX_PROPORTION);
        svgCanvas.fill(rightBottomI);
    }

    private static void drawLeftInsertIndicator(final SVGGraphics2D svgCanvas, int boxIdx)
    {
        // canvas is scaled to be in units of BASE_BOX_SIZE
        drawLeftBoxBorder(svgCanvas, boxIdx, BOUNDARY_BOX_PROPORTION);

        Rectangle2D leftTopI = new Rectangle2D.Double(boxIdx, 0.0, 1.0 / 3, BOUNDARY_BOX_PROPORTION);
        svgCanvas.fill(leftTopI);

        Rectangle2D leftBottomI = new Rectangle2D.Double(boxIdx, 1.0 - BOUNDARY_BOX_PROPORTION, 1.0 / 3, BOUNDARY_BOX_PROPORTION);
        svgCanvas.fill(leftBottomI);
    }

    private static void drawDel(SVGGraphics2D svgCanvas, int boxIdxStart, int boxIdxEnd)
    {
        // canvas is scaled to be in units of BASE_BOX_SIZE
        int delLen = boxIdxEnd - boxIdxStart + 1;

        svgCanvas.setColor(Color.BLACK);
        Rectangle2D.Double delConnector =
                new Rectangle2D.Double(boxIdxStart, 0.5 - 0.5 * DEL_CONNECTOR_BOX_PROPORTION, delLen, DEL_CONNECTOR_BOX_PROPORTION);
        svgCanvas.fill(delConnector);

        Font currentFont = svgCanvas.getFont();
        svgCanvas.setFont(INDEL_FONT);
        String delSizeStr = String.valueOf(delLen);

        // font size is based of BASE_BOX_SIZE, we do not scale font size when scaling by boxSize, so temporarily undo this scaling
        AffineTransform currentTransform = svgCanvas.getTransform();
        svgCanvas.scale(1.0 / BASE_BOX_SIZE, 1.0 / BASE_BOX_SIZE);

        double textCenterX = (boxIdxEnd + 1 - 0.5 * delLen) * BASE_BOX_SIZE;
        double textCenterY = 0.5 * BASE_BOX_SIZE;
        Rectangle2D textBackgroundRect =
                SvgUtil.scaleRectangleFromCenter(SvgUtil.getStringBoundsFromCenter(INDEL_FONT, delSizeStr, textCenterX, textCenterY), 1.5);
        svgCanvas.setColor(Color.WHITE);
        svgCanvas.fill(textBackgroundRect);

        svgCanvas.setColor(Color.BLACK);
        drawStringFromCenter(svgCanvas, delSizeStr, textCenterX, textCenterY);

        svgCanvas.setTransform(currentTransform);
        svgCanvas.setFont(currentFont);
    }

    public static SVGGraphics2D renderBaseSeq(double baseBoxSizePx, BaseRegion renderRegion, final BaseSeqViewModel bases,
            boolean shadeQuals, final Map<Integer, List<BoxBorder>> posToBordersMap, @Nullable final BaseSeqViewModel refBases)
    {
        SVGGraphics2D svgCanvas = new SVGGraphics2D(
                baseBoxSizePx * (renderRegion.baseLength() + 2 * BOX_PADDING), baseBoxSizePx);
        // work in units of BASE_BOX_SIZE
        svgCanvas.scale(baseBoxSizePx, baseBoxSizePx);

        svgCanvas.setFont(BASE_FONT);

        Integer firstBaseIdx = null;
        Integer lastBaseIdx = null;
        Integer firstOverlappingBaseIdx = null;
        Integer lastOverlappingBaseIdx = null;
        BaseViewModel prevBase = bases.getBase(renderRegion.start() - 1);
        int delLen = 0;
        for(int i = renderRegion.start(); i <= renderRegion.end(); ++i)
        {
            int boxIdx = i - renderRegion.start() + BOX_PADDING;
            BaseViewModel refBase = refBases == null ? BaseViewModel.createMissingBase() : refBases.getBase(i);
            BaseViewModel base = bases.getBase(i);

            if(!base.isMissing())
            {
                if(firstBaseIdx == null)
                {
                    firstBaseIdx = boxIdx;
                }

                lastBaseIdx = boxIdx;
            }

            if(base.isDel())
            {
                ++delLen;
            }
            else if(base.hasCharBase())
            {
                boolean matchesRef = !base.IsSoftClip && refBase.hasCharBase() && refBase.charBase() == base.charBase();

                Color bgColour = Color.WHITE;
                Color fgColour = Color.BLACK;
                if(matchesRef)
                {
                    bgColour = Color.LIGHT_GRAY;
                }
                else if(BASE_BG_COLOUR.containsKey(base.charBase()))
                {
                    bgColour = BASE_BG_COLOUR.get(base.charBase());
                    fgColour = BASE_FG_COLOUR.get(base.charBase());
                }

                if(shadeQuals)
                {
                    int baseQ = base.baseQ();
                    double lightenFactor = 1.0 - 1.0 * baseQ / MAX_BASEQ_SHADING_CUTTOFF;
                    fgColour = lighten(fgColour, lightenFactor);
                    bgColour = lighten(bgColour, lightenFactor);
                }

                svgCanvas.setColor(bgColour);
                Rectangle2D.Double box = new Rectangle2D.Double(boxIdx, 0.0, 1.0, 1.0);
                svgCanvas.fill(box);

                if(!matchesRef)
                {
                    svgCanvas.setColor(fgColour);
                    // font size is based of BASE_BOX_SIZE, we do not scale font size when scaling by boxSize, so temporarily undo this
                    // scaling
                    AffineTransform currentTransform = svgCanvas.getTransform();
                    svgCanvas.scale(1.0 / BASE_BOX_SIZE, 1.0 / BASE_BOX_SIZE);
                    drawStringFromCenter(svgCanvas, String.valueOf(base.charBase()),
                            (boxIdx + 0.5) * BASE_BOX_SIZE, 0.5 * BASE_BOX_SIZE);
                    svgCanvas.setTransform(currentTransform);
                }
            }

            // overlapping bases
            if(base.IsOverlapped)
            {
                if(firstOverlappingBaseIdx == null)
                {
                    firstOverlappingBaseIdx = boxIdx;
                }

                lastOverlappingBaseIdx = boxIdx;

                svgCanvas.setColor(OVERLAPPING_FRAGMENT_BORDER_COLOR);
                drawTopBoxBorder(svgCanvas, boxIdx, BOUNDARY_BOX_PROPORTION);
                drawBottomBoxBorder(svgCanvas, boxIdx, BOUNDARY_BOX_PROPORTION);
            }

            // del connector
            if(!base.isDel() && delLen > 0)
            {
                drawDel(svgCanvas, boxIdx - delLen, boxIdx - 1);
                delLen = 0;
            }

            // borders
            List<BoxBorder> borders = posToBordersMap.get(i);
            if(borders != null)
            {
                for(BoxBorder border : borders)
                {
                    svgCanvas.setColor(border.Col);
                    switch(border.Location)
                    {
                        case TOP:
                            drawTopBoxBorder(svgCanvas, boxIdx, 0.5 * BOUNDARY_BOX_PROPORTION);
                            break;
                        case RIGHT:
                            drawRightBoxBorder(svgCanvas, boxIdx, 0.5 * BOUNDARY_BOX_PROPORTION);
                            break;
                        case BOTTOM:
                            drawBottomBoxBorder(svgCanvas, boxIdx, 0.5 * BOUNDARY_BOX_PROPORTION);
                            break;
                        case LEFT:
                            drawLeftBoxBorder(svgCanvas, boxIdx, 0.5 * BOUNDARY_BOX_PROPORTION);
                    }
                }
            }

            // display inserts
            svgCanvas.setColor(INSERT_COLOR);
            if(base.rightInsertCount() > 0)
            {
                drawRightInsertIndicator(svgCanvas, boxIdx);
            }

            if(prevBase.rightInsertCount() > 0)
            {
                drawLeftInsertIndicator(svgCanvas, boxIdx);

                Font currentFont = svgCanvas.getFont();
                svgCanvas.setFont(INDEL_FONT);
                String insertSizeStr = String.valueOf(prevBase.rightInsertCount());

                // font size is based of BASE_BOX_SIZE, we do not scale font size when we scaled by BASE_BOX_SIZE, so temporarily undo this
                // scaling
                AffineTransform currentTransform = svgCanvas.getTransform();
                svgCanvas.scale(1.0 / BASE_BOX_SIZE, 1.0 / BASE_BOX_SIZE);

                double textCenterX = boxIdx * BASE_BOX_SIZE;
                double textCenterY = 0.5 * BASE_BOX_SIZE;
                Rectangle2D insertSizeStrRect = SvgUtil.getStringBoundsFromCenter(INDEL_FONT, insertSizeStr, textCenterX, textCenterY);
                svgCanvas.fill(insertSizeStrRect);
                svgCanvas.setColor(Color.WHITE);
                drawStringFromCenter(svgCanvas, insertSizeStr, textCenterX, textCenterY);

                svgCanvas.setTransform(currentTransform);
                svgCanvas.setFont(currentFont);
            }

            prevBase = base;
        }

        // final del connector
        if(delLen > 0)
        {
            drawDel(svgCanvas, renderRegion.end() + 1 - delLen, renderRegion.end());
        }

        // overlapping base side borders
        if(firstOverlappingBaseIdx != null && lastOverlappingBaseIdx != null)
        {
            svgCanvas.setColor(OVERLAPPING_FRAGMENT_BORDER_COLOR);
            drawLeftBoxBorder(svgCanvas, firstOverlappingBaseIdx, BOUNDARY_BOX_PROPORTION);
            drawRightBoxBorder(svgCanvas, lastOverlappingBaseIdx, BOUNDARY_BOX_PROPORTION);
        }

        if(!bases.hasOrientation() || firstBaseIdx == null || lastBaseIdx == null)
        {
            return svgCanvas;
        }

        // left orientation indicator
        if(bases.LeftIsForwardStrand)
        {
            svgCanvas.setColor(FORWARD_STRAND_COLOR);
            drawRightBoxBorder(svgCanvas, firstBaseIdx - 1, 0.5);
        }
        else
        {
            svgCanvas.setColor(REVERSE_STRAND_COLOR);
            drawReverseArrow(svgCanvas, firstBaseIdx - 0.5, 0.0, 0.5, 1.0);
        }

        // right orientation colour
        if(bases.RightIsForwardStrand)
        {
            svgCanvas.setColor(FORWARD_STRAND_COLOR);
            drawForwardArrow(svgCanvas, lastBaseIdx + 1, 0.0, 0.5, 1.0);
        }
        else
        {
            svgCanvas.setColor(REVERSE_STRAND_COLOR);
            drawLeftBoxBorder(svgCanvas, lastBaseIdx + 1, 0.5);
        }

        return svgCanvas;
    }

    public static SVGGraphics2D renderCoords(double baseBoxSizePx, final BaseRegion renderRegion, int centerPosition,
            int displayEveryNthCoord)
    {
        double scalingFactor = baseBoxSizePx / BASE_BOX_SIZE;
        double charWidth = getStringBounds(COORD_FONT, "9").getWidth();
        List<Rectangle2D> stringBounds = Lists.newArrayList();
        for(int i = renderRegion.start(); i <= renderRegion.end(); ++i)
        {
            stringBounds.add(getStringBounds(COORD_FONT, format(Locale.US, "%,d", i)));
        }

        double maxStringWidth = stringBounds.stream().mapToDouble(Rectangle2D::getWidth).max().getAsDouble();

        SVGGraphics2D svgCanvas = new SVGGraphics2D(
                scalingFactor * BASE_BOX_SIZE * (renderRegion.baseLength() + 2 * BOX_PADDING),
                scalingFactor * (maxStringWidth + 2 * charWidth));
        svgCanvas.setFont(COORD_FONT);

        Rectangle2D.Double backgroundRect = new Rectangle2D.Double(0.0, 0.0, svgCanvas.getWidth(), svgCanvas.getHeight());
        svgCanvas.setColor(Color.LIGHT_GRAY);
        svgCanvas.fill(backgroundRect);

        svgCanvas.scale(scalingFactor, scalingFactor);
        AffineTransform scalingTransform = svgCanvas.getTransform();

        svgCanvas.setColor(Color.BLACK);
        for(int i = renderRegion.start(); i <= renderRegion.end(); ++i)
        {
            svgCanvas.setTransform(scalingTransform);

            int boxIdx = i - renderRegion.start() + BOX_PADDING;
            int basePosOffset = i - centerPosition;
            boolean displayCoord = basePosOffset % displayEveryNthCoord == 0;
            if(!displayCoord)
            {
                Rectangle2D smallTick = new Rectangle2D.Double(
                        BASE_BOX_SIZE * (boxIdx + 0.5 - 0.25 * BOUNDARY_BOX_PROPORTION),
                        svgCanvas.getHeight() / scalingFactor - 0.75 * charWidth,
                        BASE_BOX_SIZE * 0.5 * BOUNDARY_BOX_PROPORTION,
                        0.75 * charWidth);
                svgCanvas.fill(smallTick);
                continue;
            }

            Rectangle2D largeTick = new Rectangle2D.Double(
                    BASE_BOX_SIZE * (boxIdx + 0.5 - 0.25 * BOUNDARY_BOX_PROPORTION),
                    svgCanvas.getHeight() / scalingFactor - 1.5 * charWidth,
                    BASE_BOX_SIZE * 0.5 * BOUNDARY_BOX_PROPORTION,
                    1.5 * charWidth);
            svgCanvas.fill(largeTick);

            Rectangle2D stringBound = stringBounds.get(i - renderRegion.start());

            svgCanvas.rotate(Math.PI / 2);
            svgCanvas.translate(maxStringWidth - 0.5 * stringBound.getWidth(), -BASE_BOX_SIZE * (boxIdx + 0.5));
            drawStringFromCenter(svgCanvas, format(Locale.US, "%,d", i), 0, 0);
        }

        return svgCanvas;
    }

    public static void drawForwardArrow(final SVGGraphics2D svgCanvas, double left, double top, double width, double height)
    {
        Path2D.Double forwardArrowPath = new Path2D.Double();
        forwardArrowPath.moveTo(left, top);
        forwardArrowPath.lineTo(left + width, top + 0.5 * height);
        forwardArrowPath.lineTo(left, top + height);
        forwardArrowPath.closePath();

        svgCanvas.fill(forwardArrowPath);
    }

    public static void drawReverseArrow(final SVGGraphics2D svgCanvas, double left, double top, double width, double height)
    {
        Path2D.Double reverseArrowPath = new Path2D.Double();
        reverseArrowPath.moveTo(left, top + 0.5 * height);
        reverseArrowPath.lineTo(left + width, top);
        reverseArrowPath.lineTo(left + width, top + height);
        reverseArrowPath.closePath();

        svgCanvas.fill(reverseArrowPath);
    }

    public static SVGGraphics2D renderColoredBox(double sizePx, final Color color)
    {
        SVGGraphics2D svgCanvas = new SVGGraphics2D(sizePx, sizePx);
        svgCanvas.setColor(color);
        Rectangle2D rect = new Rectangle2D.Double(0.0, 0.0, sizePx, sizePx);
        svgCanvas.fill(rect);
        return svgCanvas;
    }

    public enum BorderLocation
    {
        TOP,
        RIGHT,
        BOTTOM,
        LEFT
    }

    public static class BoxBorder
    {
        public final BorderLocation Location;
        public final Color Col;

        public BoxBorder(final BorderLocation location, final Color col)
        {
            Location = location;
            Col = col;
        }
    }
}
