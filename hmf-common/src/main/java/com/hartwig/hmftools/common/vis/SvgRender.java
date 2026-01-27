package com.hartwig.hmftools.common.vis;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;
import static java.util.Map.entry;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions.enforceChrPrefix;
import static com.hartwig.hmftools.common.vis.ColorUtil.PURPLE;
import static com.hartwig.hmftools.common.vis.ColorUtil.lighten;
import static com.hartwig.hmftools.common.vis.SvgUtil.Alignment.CENTER;
import static com.hartwig.hmftools.common.vis.SvgUtil.Alignment.LEFT;
import static com.hartwig.hmftools.common.vis.SvgUtil.Alignment.RIGHT;
import static com.hartwig.hmftools.common.vis.SvgUtil.drawStringFromCenter;
import static com.hartwig.hmftools.common.vis.SvgUtil.getStringBounds;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.Shape;
import java.awt.Stroke;
import java.awt.geom.AffineTransform;
import java.awt.geom.Line2D;
import java.awt.geom.Path2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.util.List;
import java.util.Locale;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.vis.SvgUtil.Alignment;

import org.jetbrains.annotations.Nullable;
import org.jfree.svg.SVGGraphics2D;

public final class SvgRender
{
    public static final Color FORWARD_STRAND_COLOR = new Color(81, 144, 207);
    public static final Color REVERSE_STRAND_COLOR = new Color(242, 167, 121);
    public static final Color OVERLAPPING_FRAGMENT_BORDER_COLOR = lighten(PURPLE, 0.5f);

    private static final Point2D.Double ZERO_2D = new Point2D.Double(0.0, 0.0);

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

    private static final Color AMINO_ACID_MATCH_LIGHT_BG = new Color(92, 92, 164);
    private static final Color AMINO_ACID_MATCH_DARK_BG = new Color(12, 12, 120);
    private static final List<Color> AMINO_ACID_MATCH_BG_COLOURS = Lists.newArrayList(AMINO_ACID_MATCH_LIGHT_BG, AMINO_ACID_MATCH_DARK_BG);

    private static final Color AMINO_ACID_MISMATCH_LIGHT_BG = new Color(179, 130, 77);
    private static final Color AMINO_ACID_MISMATCH_DARK_BG = new Color(182, 103, 18);
    private static final List<Color> AMINO_ACID_MISMATCH_BG_COLOURS = Lists.newArrayList(
            AMINO_ACID_MISMATCH_LIGHT_BG, AMINO_ACID_MISMATCH_DARK_BG);

    public static final int BOX_PADDING = 1;
    public static final double BASE_BOX_SIZE = 30.0;
    private static final double DEL_CONNECTOR_BOX_PROPORTION = 0.2;
    private static final double BOUNDARY_BOX_PROPORTION = 0.2;

    // font sizes are tied to BASE_BOX_SIZE
    private static final int INDEL_FONT_SIZE = 18;
    private static final int BASE_FONT_SIZE = 26;
    private static final int COORD_FONT_SIZE = 30;
    private static final Font INDEL_FONT = new Font("SansSerif", Font.PLAIN, INDEL_FONT_SIZE);
    private static final Font BASE_FONT = new Font("SansSerif", Font.PLAIN, BASE_FONT_SIZE);
    public static final Font COORD_FONT = new Font("SansSerif", Font.PLAIN, COORD_FONT_SIZE);

    private static final int MAX_BASEQ_SHADING_CUTTOFF = 37;

    private static final String INSERT_LABEL = "INS";

    private SvgRender() {}

    private static void drawTopBoxBorder(final SVGGraphics2D svgCanvas, final Point2D.Double boxOffset, int boxIdx, double boxPropWidth)
    {
        // canvas is scaled to be in units of BASE_BOX_SIZE
        Rectangle2D rect = new Rectangle2D.Double(boxIdx + boxOffset.x, boxOffset.y, 1.0, boxPropWidth);
        svgCanvas.fill(rect);
    }

    private static void drawRightBoxBorder(final SVGGraphics2D svgCanvas, final Point2D.Double boxOffset, int boxIdx, double boxPropWidth)
    {
        // canvas is scaled to be in units of BASE_BOX_SIZE
        Rectangle2D rect = new Rectangle2D.Double(boxIdx + 1 - boxPropWidth + boxOffset.x, boxOffset.y, boxPropWidth, 1.0);
        svgCanvas.fill(rect);
    }

    private static void drawBottomBoxBorder(final SVGGraphics2D svgCanvas, final Point2D.Double boxOffset, int boxIdx, double boxPropWidth)
    {
        // canvas is scaled to be in units of BASE_BOX_SIZE
        Rectangle2D rect = new Rectangle2D.Double(boxIdx + boxOffset.x, 1.0 - boxPropWidth + boxOffset.y, 1.0, boxPropWidth);
        svgCanvas.fill(rect);
    }

    private static void drawLeftBoxBorder(final SVGGraphics2D svgCanvas, final Point2D.Double boxOffset, int boxIdx, double boxPropWidth)
    {
        // canvas is scaled to be in units of BASE_BOX_SIZE
        Rectangle2D rect = new Rectangle2D.Double(boxIdx + boxOffset.x, boxOffset.y, boxPropWidth, 1.0);
        svgCanvas.fill(rect);
    }

    private static void drawRightInsertIndicator(final SVGGraphics2D svgCanvas, final Point2D.Double boxOffset, int boxIdx)
    {
        // canvas is scaled to be in units of BASE_BOX_SIZE
        drawRightBoxBorder(svgCanvas, boxOffset, boxIdx, BOUNDARY_BOX_PROPORTION);

        Rectangle2D rightTopI = new Rectangle2D.Double(boxIdx + 1 - 1.0 / 3 + boxOffset.x, boxOffset.y, 1.0 / 3, BOUNDARY_BOX_PROPORTION);
        svgCanvas.fill(rightTopI);

        Rectangle2D rightBottomI = new Rectangle2D.Double(
                boxIdx + 1 - 1.0 / 3 + boxOffset.x, 1.0 - BOUNDARY_BOX_PROPORTION + boxOffset.y, 1.0 / 3, BOUNDARY_BOX_PROPORTION);
        svgCanvas.fill(rightBottomI);
    }

    private static void drawLeftInsertIndicator(final SVGGraphics2D svgCanvas, final Point2D.Double boxOffset, int boxIdx)
    {
        // canvas is scaled to be in units of BASE_BOX_SIZE
        drawLeftBoxBorder(svgCanvas, boxOffset, boxIdx, BOUNDARY_BOX_PROPORTION);

        Rectangle2D leftTopI = new Rectangle2D.Double(boxIdx + boxOffset.x, boxOffset.y, 1.0 / 3, BOUNDARY_BOX_PROPORTION);
        svgCanvas.fill(leftTopI);

        Rectangle2D leftBottomI = new Rectangle2D.Double(
                boxIdx + boxOffset.x, 1.0 - BOUNDARY_BOX_PROPORTION + boxOffset.y, 1.0 / 3, BOUNDARY_BOX_PROPORTION);
        svgCanvas.fill(leftBottomI);
    }

    private static void drawDel(final SVGGraphics2D svgCanvas, final Point2D.Double boxOffset, @Nullable final Integer labelDelLen,
            int boxIdxStart, int boxIdxEnd)
    {
        // canvas is scaled to be in units of BASE_BOX_SIZE
        int delLen = boxIdxEnd - boxIdxStart + 1;

        svgCanvas.setColor(Color.BLACK);
        Rectangle2D.Double delConnector =
                new Rectangle2D.Double(
                        boxIdxStart + boxOffset.x,
                        0.5 - 0.5 * DEL_CONNECTOR_BOX_PROPORTION + boxOffset.y, delLen, DEL_CONNECTOR_BOX_PROPORTION);
        svgCanvas.fill(delConnector);

        if(labelDelLen == null)
            return;

        String delSizeStr = String.valueOf(labelDelLen);
        renderText(svgCanvas, boxOffset, boxIdxEnd + 1 - 0.5 * delLen, 0.5, INDEL_FONT, Color.BLACK, delSizeStr, Color.WHITE, 1.5, CENTER);
    }

    private static void renderText(final SVGGraphics2D svgCanvas, final Point2D.Double boxOffset, double boxX, double boxY,
            final Font font, final Color fgColor, final String text, @Nullable final Color bgColor, double backgroundBoxExpandFactor, final Alignment alignment)
    {
        Font currentFont = svgCanvas.getFont();
        svgCanvas.setFont(font);

        // font size is based of BASE_BOX_SIZE, we do not scale font size when scaling by boxSize, so temporarily undo this scaling
        AffineTransform currentTransform = svgCanvas.getTransform();
        svgCanvas.scale(1.0 / BASE_BOX_SIZE, 1.0 / BASE_BOX_SIZE);

        Rectangle2D stringBounds = null;
        if(alignment == CENTER)
        {

            double textCenterX = (boxX + boxOffset.x) * BASE_BOX_SIZE;
            double textCenterY = (boxY + boxOffset.y) * BASE_BOX_SIZE;

            if(bgColor != null)
            {
                Rectangle2D textBackgroundRect = SvgUtil.scaleRectangleFromCenter(
                        SvgUtil.getStringBoundsFromCenter(font, text, textCenterX, textCenterY), backgroundBoxExpandFactor);
                svgCanvas.setColor(bgColor);
                svgCanvas.fill(textBackgroundRect);
            }

            svgCanvas.setColor(fgColor);
            drawStringFromCenter(svgCanvas, text, textCenterX, textCenterY);
        }
        else
        {
            stringBounds = getStringBounds(font, text);
        }

        svgCanvas.setTransform(currentTransform);
        svgCanvas.setFont(currentFont);

        if(alignment == LEFT)
        {
            double boxCenterX = boxX + stringBounds.getCenterX() / BASE_BOX_SIZE;
            renderText(svgCanvas, boxOffset, boxCenterX, boxY, font, fgColor, text, bgColor, backgroundBoxExpandFactor, CENTER);
        }
        else if(alignment == RIGHT)
        {
            double boxCenterX = boxX - stringBounds.getCenterX() / BASE_BOX_SIZE;
            renderText(svgCanvas, boxOffset, boxCenterX, boxY, font, fgColor, text, bgColor, backgroundBoxExpandFactor, CENTER);
        }
    }

    public static SVGGraphics2D renderBaseSeq(double baseBoxSizePx, final BaseRegion renderRegion, final BaseSeqViewModel bases,
            boolean shadeQuals, final Map<Integer, List<BoxBorder>> posToBordersMap, @Nullable final BaseSeqViewModel refBases)
    {
        SVGGraphics2D svgCanvas = new SVGGraphics2D(
                baseBoxSizePx * (renderRegion.baseLength() + 2 * BOX_PADDING), baseBoxSizePx);
        renderBaseSeq(svgCanvas, ZERO_2D, baseBoxSizePx, renderRegion, bases, shadeQuals, posToBordersMap, refBases, true, true, false);
        return svgCanvas;
    }

    public static void renderBaseSeq(final SVGGraphics2D svgCanvas, final Point2D.Double boxOffset, double baseBoxSizePx,
            final BaseRegion renderRegion, final BaseSeqViewModel bases, boolean shadeQuals,
            final Map<Integer, List<BoxBorder>> posToBordersMap, @Nullable final BaseSeqViewModel refBases)
    {
        renderBaseSeq(svgCanvas, boxOffset, baseBoxSizePx, renderRegion, bases, shadeQuals, posToBordersMap, refBases, true, true, false);
    }

    public static void renderBaseSeq(final SVGGraphics2D svgCanvas, final Point2D.Double boxOffset, double baseBoxSizePx,
            final BaseRegion renderRegion, final BaseSeqViewModel bases, boolean shadeQuals,
            final Map<Integer, List<BoxBorder>> posToBordersMap, @Nullable final BaseSeqViewModel refBases,
            boolean renderLeftOrientationMarker, boolean renderRightOrientationMarker, boolean reverseOrientationMarkers)
    {
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
                boolean matchesRef = !base.isSoftClip() && refBase.hasCharBase() && refBase.charBase() == base.charBase();

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
                Rectangle2D.Double box = new Rectangle2D.Double(boxIdx + boxOffset.x, boxOffset.y, 1.0, 1.0);
                svgCanvas.fill(box);

                if(!matchesRef)
                    renderText(svgCanvas, boxOffset,
                            boxIdx + 0.5, 0.5, svgCanvas.getFont(), fgColour, String.valueOf(base.charBase()), null, 1.0, CENTER);
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
                drawTopBoxBorder(svgCanvas, boxOffset, boxIdx, BOUNDARY_BOX_PROPORTION);
                drawBottomBoxBorder(svgCanvas, boxOffset, boxIdx, BOUNDARY_BOX_PROPORTION);
            }

            // del connector
            if(!base.isDel() && delLen > 0)
            {
                drawDel(svgCanvas, boxOffset, delLen, boxIdx - delLen, boxIdx - 1);
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
                            drawTopBoxBorder(svgCanvas, boxOffset, boxIdx, 0.5 * BOUNDARY_BOX_PROPORTION);
                            break;
                        case RIGHT:
                            drawRightBoxBorder(svgCanvas, boxOffset, boxIdx, 0.5 * BOUNDARY_BOX_PROPORTION);
                            break;
                        case BOTTOM:
                            drawBottomBoxBorder(svgCanvas, boxOffset, boxIdx, 0.5 * BOUNDARY_BOX_PROPORTION);
                            break;
                        case LEFT:
                            drawLeftBoxBorder(svgCanvas, boxOffset, boxIdx, 0.5 * BOUNDARY_BOX_PROPORTION);
                    }
                }
            }

            // display inserts
            svgCanvas.setColor(INSERT_COLOR);
            if(base.rightInsertCount() > 0)
            {
                drawRightInsertIndicator(svgCanvas, boxOffset, boxIdx);
            }

            if(prevBase.rightInsertCount() > 0)
            {
                drawLeftInsertIndicator(svgCanvas, boxOffset, boxIdx);
                String insertSizeStr = String.valueOf(prevBase.rightInsertCount());
                renderText(svgCanvas, boxOffset, boxIdx, 0.5, INDEL_FONT, Color.WHITE, insertSizeStr, INSERT_COLOR, 1.0, CENTER);
            }

            prevBase = base;
        }

        // final del connector
        if(delLen > 0)
        {
            int boxEndIdx = renderRegion.end() - renderRegion.start() + BOX_PADDING;
            drawDel(svgCanvas, boxOffset, delLen, boxEndIdx + 1 - delLen, boxEndIdx);
        }

        // overlapping base side borders
        if(firstOverlappingBaseIdx != null && lastOverlappingBaseIdx != null)
        {
            svgCanvas.setColor(OVERLAPPING_FRAGMENT_BORDER_COLOR);
            drawLeftBoxBorder(svgCanvas, boxOffset, firstOverlappingBaseIdx, BOUNDARY_BOX_PROPORTION);
            drawRightBoxBorder(svgCanvas, boxOffset, lastOverlappingBaseIdx, BOUNDARY_BOX_PROPORTION);
        }

        if(!bases.hasOrientation() || firstBaseIdx == null || lastBaseIdx == null)
            return;

        // left orientation indicator
        if(renderLeftOrientationMarker)
        {
            boolean leftIsForwardStrand = bases.LeftIsForwardStrand;
            if(reverseOrientationMarkers)
                leftIsForwardStrand = !leftIsForwardStrand;

            if(leftIsForwardStrand)
            {
                svgCanvas.setColor(FORWARD_STRAND_COLOR);
                drawRightBoxBorder(svgCanvas, boxOffset, firstBaseIdx - 1, 0.5);
            }
            else
            {
                svgCanvas.setColor(REVERSE_STRAND_COLOR);
                drawReverseArrow(svgCanvas, boxOffset, firstBaseIdx - 0.5, 0.0, 0.5, 1.0);
            }
        }

        // right orientation colour
        if(renderRightOrientationMarker)
        {
            boolean rightIsForwardStrand = bases.RightIsForwardStrand;
            if(reverseOrientationMarkers)
                rightIsForwardStrand = !rightIsForwardStrand;

            if(rightIsForwardStrand)
            {
                svgCanvas.setColor(FORWARD_STRAND_COLOR);
                drawForwardArrow(svgCanvas, boxOffset, lastBaseIdx + 1, 0.0, 0.5, 1.0);
            }
            else
            {
                svgCanvas.setColor(REVERSE_STRAND_COLOR);
                drawLeftBoxBorder(svgCanvas, boxOffset, lastBaseIdx + 1, 0.5);
            }
        }
    }

    public static SVGGraphics2D renderCoords(double baseBoxSizePx, final BaseRegion renderRegion, int centerPosition,
            int displayEveryNthCoord)
    {
        return renderCoords(baseBoxSizePx, renderRegion, centerPosition, displayEveryNthCoord, false);
    }

    public static SVGGraphics2D renderCoords(double baseBoxSizePx, final BaseRegion renderRegion, int centerPosition,
            int displayEveryNthCoord, boolean reverse)
    {
        double scalingFactor = baseBoxSizePx / BASE_BOX_SIZE;
        double charWidth = getStringBounds(COORD_FONT, "9").getWidth();
        List<Rectangle2D> stringBounds = Lists.newArrayList();
        for(int i = renderRegion.start(); i <= renderRegion.end(); ++i)
            stringBounds.add(getStringBounds(COORD_FONT, format(Locale.US, "%,d", i)));

        double maxStringWidth = stringBounds.stream().mapToDouble(Rectangle2D::getWidth).max().getAsDouble();

        SVGGraphics2D svgCanvas = new SVGGraphics2D(
                scalingFactor * BASE_BOX_SIZE * (renderRegion.baseLength() + 2 * BOX_PADDING),
                scalingFactor * (maxStringWidth + 2 * charWidth));
        Point2D.Double canvasSize = new Point2D.Double(svgCanvas.getWidth(), svgCanvas.getHeight());

        renderCoords(svgCanvas, ZERO_2D, canvasSize, baseBoxSizePx, renderRegion, centerPosition, displayEveryNthCoord, maxStringWidth, reverse);
        return svgCanvas;
    }

    public static void renderCoords(final SVGGraphics2D svgCanvas, final Point2D.Double boxOffset, final Point2D.Double canvasSize,
            double baseBoxSizePx, final BaseRegion renderRegion, int centerPosition, int displayEveryNthCoord, double maxStringWidth, boolean reverse)
    {
        double scalingFactor = baseBoxSizePx / BASE_BOX_SIZE;
        double charWidth = getStringBounds(COORD_FONT, "9").getWidth();
        List<Rectangle2D> stringBounds = Lists.newArrayList();
        for(int i = renderRegion.start(); i <= renderRegion.end(); ++i)
            stringBounds.add(getStringBounds(COORD_FONT, format(Locale.US, "%,d", i)));

        svgCanvas.setFont(COORD_FONT);

        Rectangle2D.Double backgroundRect = new Rectangle2D.Double(
                baseBoxSizePx * boxOffset.x,
                baseBoxSizePx * boxOffset.y, canvasSize.x + baseBoxSizePx * boxOffset.x, canvasSize.y + baseBoxSizePx * boxOffset.y);
        svgCanvas.setColor(Color.LIGHT_GRAY);
        svgCanvas.fill(backgroundRect);

        svgCanvas.scale(scalingFactor, scalingFactor);
        AffineTransform scalingTransform = svgCanvas.getTransform();

        svgCanvas.setColor(Color.BLACK);
        for(int box = 0; box < renderRegion.baseLength(); ++box)
        {
            svgCanvas.setTransform(scalingTransform);

            int boxIdx = box + BOX_PADDING;
            int basePos = reverse ? renderRegion.end() - box : renderRegion.start() + box;
            int basePosOffset = basePos - centerPosition;
            boolean displayCoord = basePosOffset % displayEveryNthCoord == 0;
            if(!displayCoord)
            {
                Rectangle2D smallTick = new Rectangle2D.Double(
                        BASE_BOX_SIZE * (boxIdx + 0.5 - 0.25 * BOUNDARY_BOX_PROPORTION + boxOffset.x),
                        canvasSize.y / scalingFactor - 0.75 * charWidth + BASE_BOX_SIZE * boxOffset.y,
                        BASE_BOX_SIZE * 0.5 * BOUNDARY_BOX_PROPORTION,
                        0.75 * charWidth);
                svgCanvas.fill(smallTick);
                continue;
            }

            Rectangle2D largeTick = new Rectangle2D.Double(
                    BASE_BOX_SIZE * (boxIdx + 0.5 - 0.25 * BOUNDARY_BOX_PROPORTION + boxOffset.x),
                    canvasSize.y / scalingFactor - 1.5 * charWidth + BASE_BOX_SIZE * boxOffset.y,
                    BASE_BOX_SIZE * 0.5 * BOUNDARY_BOX_PROPORTION,
                    1.5 * charWidth);
            svgCanvas.fill(largeTick);

            int stringIdx = reverse ? stringBounds.size() - 1 - box : box;
            Rectangle2D stringBound = stringBounds.get(stringIdx);

            svgCanvas.rotate(Math.PI / 2);
            svgCanvas.translate(
                    maxStringWidth - 0.5 * stringBound.getWidth() + BASE_BOX_SIZE * boxOffset.y,
                    -BASE_BOX_SIZE * (boxIdx + 0.5 + boxOffset.x));
            drawStringFromCenter(svgCanvas, format(Locale.US, "%,d", basePos), 0, 0);
        }
    }

    public record ChrLabel(String chromosome, int position, BaseRegion viewRegion, boolean isReversed, Alignment alignment) {}

    public static SVGGraphics2D renderChrLabels(double baseBoxSizePx, final List<ChrLabel> labels)
    {
        int startIdx = labels.get(0).viewRegion.start();
        int endIdx = labels.get(labels.size() - 1).viewRegion.end();
        int boxLength = endIdx - startIdx + 1;
        SVGGraphics2D svgCanvas = new SVGGraphics2D(baseBoxSizePx * boxLength, baseBoxSizePx);

        // work in units of BASE_BOX_SIZE
        svgCanvas.scale(baseBoxSizePx, baseBoxSizePx);

        svgCanvas.setFont(BASE_FONT);
        for(ChrLabel label : labels)
        {
            String chromosome = label.chromosome;
            BaseRegion viewRegion = label.viewRegion;

            double boxX;
            String fullLabel;
            if(chromosome == null)
            {
                fullLabel = INSERT_LABEL;
                boxX = viewRegion.start() + 0.5d * viewRegion.baseLength();
                renderText(svgCanvas, ZERO_2D, boxX, 0.5, svgCanvas.getFont(), Color.BLACK, fullLabel, null, 1.0, CENTER);
                continue;
            }

            chromosome = enforceChrPrefix(chromosome);
            Alignment alignment = label.alignment;
            String strandLabel = label.isReversed ? "<<" : ">>";
            String chrPositionLabel = chromosome + ":" + label.position;
            if(alignment == CENTER)
            {
                fullLabel = strandLabel + " " + chrPositionLabel + " " + strandLabel;
                boxX = viewRegion.start() + 0.5d * viewRegion.baseLength();
            }
            else if(alignment == LEFT)
            {
                fullLabel = chrPositionLabel + " " + strandLabel;
                boxX = viewRegion.start() + BOX_PADDING;
            }
            else
            {
                fullLabel = strandLabel + " " + chrPositionLabel;
                boxX = viewRegion.start() + viewRegion.baseLength();
            }

            renderText(svgCanvas, ZERO_2D, boxX, 0.5, svgCanvas.getFont(), Color.BLACK, fullLabel, null, 1.0, alignment);
        }

        return svgCanvas;
    }

    private static void drawForwardArrow(final SVGGraphics2D svgCanvas, final Point2D.Double offset, double left, double top, double width,
            double height)
    {
        Path2D.Double forwardArrowPath = new Path2D.Double();
        forwardArrowPath.moveTo(left + offset.x, top + offset.y);
        forwardArrowPath.lineTo(left + width + offset.x, top + 0.5 * height + offset.y);
        forwardArrowPath.lineTo(left + offset.x, top + height + offset.y);
        forwardArrowPath.closePath();

        svgCanvas.fill(forwardArrowPath);
    }

    private static void drawReverseArrow(final SVGGraphics2D svgCanvas, final Point2D.Double offset, double left, double top, double width,
            double height)
    {
        Path2D.Double reverseArrowPath = new Path2D.Double();
        reverseArrowPath.moveTo(left + offset.x, top + 0.5 * height + offset.y);
        reverseArrowPath.lineTo(left + width + offset.x, top + offset.y);
        reverseArrowPath.lineTo(left + width + offset.x, top + height + offset.y);
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

    public record RenderedGeneData(SVGGraphics2D refSvgCanvas, SVGGraphics2D altSvgCanvas) {}

    public static RenderedGeneData renderGeneData(double baseBoxSizePx, final BaseRegion renderRegion, boolean posStrand,
            final List<GeneRegionViewModel> refViewModels, final List<GeneRegionViewModel> altViewModels)
    {
        SVGGraphics2D refSvgCanvas = new SVGGraphics2D(
                baseBoxSizePx * (renderRegion.baseLength() + 2 * BOX_PADDING), baseBoxSizePx);
        renderGeneRegions(refSvgCanvas, baseBoxSizePx, renderRegion, posStrand, refViewModels, true);

        SVGGraphics2D altSvgCanvas = new SVGGraphics2D(
                baseBoxSizePx * (renderRegion.baseLength() + 2 * BOX_PADDING), baseBoxSizePx);
        renderGeneRegions(altSvgCanvas, baseBoxSizePx, renderRegion, posStrand, altViewModels, false);

        return new RenderedGeneData(refSvgCanvas, altSvgCanvas);
    }

    private static void renderGeneRegions(final SVGGraphics2D svgCanvas, double baseBoxSizePx, final BaseRegion renderRegion,
            boolean posStrand, final List<GeneRegionViewModel> viewModels, boolean renderRef)
    {
        // work in units of BASE_BOX_SIZE
        svgCanvas.scale(baseBoxSizePx, baseBoxSizePx);
        svgCanvas.setFont(BASE_FONT);
        for(int i = 0; i < viewModels.size(); i++)
        {
            GeneRegionViewModel viewModel = viewModels.get(i);
            if(viewModel.region().start() > renderRegion.end())
                break;

            if(viewModel.region().end() < renderRegion.start())
                continue;

            if(viewModel instanceof GeneRegionViewModel.IntronicRegionViewModel)
            {
                svgCanvas.setColor(AMINO_ACID_MATCH_LIGHT_BG);
                Stroke currentStroke = svgCanvas.getStroke();
                Stroke stroke = new BasicStroke((float) (2.0 / BASE_BOX_SIZE));
                svgCanvas.setStroke(stroke);
                double left = max(viewModel.region().start(), renderRegion.start()) - renderRegion.start() + BOX_PADDING;
                double right = min(viewModel.region().end(), renderRegion.end()) + 1 - renderRegion.start() + BOX_PADDING;
                Shape line = new Line2D.Double(left, 0.5, right, 0.5);
                svgCanvas.draw(line);
                svgCanvas.setStroke(currentStroke);
                continue;
            }

            if(viewModel instanceof GeneRegionViewModel.NonCodingExonicRegionViewModel)
            {
                svgCanvas.setColor(AMINO_ACID_MATCH_LIGHT_BG);
                double left = max(viewModel.region().start(), renderRegion.start()) - renderRegion.start() + BOX_PADDING;
                double right = min(viewModel.region().end(), renderRegion.end()) - renderRegion.start() + BOX_PADDING;
                double len = right - left + 1;
                Rectangle2D.Double box = new Rectangle2D.Double(left, 0.25, len, 0.5);
                svgCanvas.fill(box);
                continue;
            }

            BaseRegion region = viewModel.region();
            int startPos = max(region.start(), renderRegion.start());
            int endPos = min(region.end(), renderRegion.end());
            int boxWidth = endPos - startPos + 1;
            int boxIdx = startPos - renderRegion.start() + BOX_PADDING;
            if(viewModel instanceof GeneRegionViewModel.DelViewModel)
            {
                drawDel(svgCanvas, ZERO_2D, null, boxIdx, boxIdx + boxWidth - 1);
                continue;
            }

            GeneRegionViewModel.AminoAcidViewModel aaModel = (GeneRegionViewModel.AminoAcidViewModel) viewModel;
            char aa = renderRef ? aaModel.ref() : aaModel.alt();
            boolean matchesRef = renderRef || aaModel.matchesRef();

            Color bgColour;
            if(aa == '*')
                bgColour = Color.RED;
            else if(aaModel.isStart())
                bgColour = Color.GREEN;
            else if(matchesRef)
                bgColour = AMINO_ACID_MATCH_BG_COLOURS.get(i % AMINO_ACID_MATCH_BG_COLOURS.size());
            else
                bgColour = AMINO_ACID_MISMATCH_BG_COLOURS.get(i % AMINO_ACID_MISMATCH_BG_COLOURS.size());

            svgCanvas.setColor(bgColour);
            Rectangle2D.Double box = new Rectangle2D.Double(boxIdx, 0.0, boxWidth, 1.0);
            svgCanvas.fill(box);

            Color fgColour = Color.WHITE;
            renderText(svgCanvas, ZERO_2D, boxIdx + 0.5 * boxWidth, 0.5, svgCanvas.getFont(), fgColour,
		    String.valueOf(aa), null, 1.0, CENTER);
        }

        // left orientation indicator
        svgCanvas.setColor(Color.BLUE);
        if(posStrand)
        {
            drawForwardArrow(svgCanvas, ZERO_2D, 0.5, 0.0, 0.5, 1.0);
            drawForwardArrow(svgCanvas, ZERO_2D, renderRegion.baseLength() + BOX_PADDING, 0.0, 0.5, 1.0);
        }
        else
        {
            drawReverseArrow(svgCanvas, ZERO_2D, 0.5, 0.0, 0.5, 1.0);
            drawReverseArrow(svgCanvas, ZERO_2D, renderRegion.baseLength() + BOX_PADDING, 0.0, 0.5, 1.0);
        }

        // display inserts
        if(renderRef)
            return;

        for(GeneRegionViewModel viewModel : viewModels)
        {
            if(viewModel.region().start() > renderRegion.end())
                break;

            if(viewModel.region().end() < renderRegion.start())
                continue;

            if(viewModel instanceof GeneRegionViewModel.DelViewModel)
                continue;

            if(viewModel.basesBeforeRightInsert() <= 0)
                continue;

            BaseRegion region = viewModel.region();
            int startPos = max(region.start(), renderRegion.start());
            int endPos = min(region.end(), renderRegion.end());
            int boxWidth = endPos - startPos + 1;
            int boxIdx = startPos - renderRegion.start() + BOX_PADDING;

            int indelPos = posStrand
                    ? region.start() + viewModel.basesBeforeRightInsert() - 1
                    : viewModel.region().end() - viewModel.basesBeforeRightInsert() + 1;
            if(!renderRegion.containsPosition(indelPos))
                continue;

            int indelBoxIdx = indelPos - renderRegion.start() + BOX_PADDING;
            svgCanvas.setColor(INSERT_COLOR);
            if(posStrand)
            {
                drawRightInsertIndicator(svgCanvas, ZERO_2D, indelBoxIdx);
                drawLeftInsertIndicator(svgCanvas, ZERO_2D, indelBoxIdx + 1);
            }
            else
            {
                drawLeftInsertIndicator(svgCanvas, ZERO_2D, indelBoxIdx);
                drawRightInsertIndicator(svgCanvas, ZERO_2D, indelBoxIdx - 1);
            }

            // re-render alt amino acid
            if(!(viewModel instanceof final GeneRegionViewModel.AminoAcidViewModel aaModel))
                continue;

            Color fgColour = Color.WHITE;
            renderText(svgCanvas, ZERO_2D,
                    boxIdx + 0.5 * boxWidth, 0.5, svgCanvas.getFont(), fgColour, String.valueOf(aaModel.alt()), null, 1.0, CENTER);
        }
    }
}
