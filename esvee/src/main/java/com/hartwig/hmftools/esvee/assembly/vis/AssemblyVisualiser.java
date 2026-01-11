package com.hartwig.hmftools.esvee.assembly.vis;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.vis.HtmlUtil.BASE_FONT_STYLE;
import static com.hartwig.hmftools.common.vis.HtmlUtil.JQUERY_SCRIPT;
import static com.hartwig.hmftools.common.vis.HtmlUtil.getJavascript;
import static com.hartwig.hmftools.common.vis.HtmlUtil.styledTable;
import static com.hartwig.hmftools.common.vis.SvgRender.BASE_BOX_SIZE;
import static com.hartwig.hmftools.common.vis.SvgRender.BOX_PADDING;
import static com.hartwig.hmftools.common.vis.SvgRender.COORD_FONT;
import static com.hartwig.hmftools.common.vis.SvgRender.renderBaseSeq;
import static com.hartwig.hmftools.common.vis.SvgUtil.getStringBounds;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.assembly.vis.AssemblyVisConstants.DISPLAY_EVERY_NTH_COORD;
import static com.hartwig.hmftools.esvee.assembly.vis.AssemblyVisConstants.READ_HEIGHT_PX;
import static com.hartwig.hmftools.esvee.assembly.vis.AssemblyVisConstants.VIEW_REGION_SIZE;
import static com.hartwig.hmftools.esvee.assembly.vis.AssemblyVisConstants.VIS_DIR;

import static j2html.TagCreator.body;
import static j2html.TagCreator.div;
import static j2html.TagCreator.header;
import static j2html.TagCreator.html;
import static j2html.TagCreator.rawHtml;
import static j2html.TagCreator.td;
import static j2html.TagCreator.tr;

import java.awt.Color;
import java.awt.geom.AffineTransform;
import java.awt.geom.Point2D;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayDeque;
import java.util.List;
import java.util.Locale;
import java.util.Queue;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.vis.BaseSeqViewModel;
import com.hartwig.hmftools.common.vis.CssBuilder;
import com.hartwig.hmftools.common.vis.CssSize;
import com.hartwig.hmftools.common.vis.SvgRender;
import com.hartwig.hmftools.esvee.assembly.AssemblyConfig;
import com.hartwig.hmftools.esvee.assembly.alignment.AlignData;
import com.hartwig.hmftools.esvee.assembly.alignment.AssemblyAlignment;
import com.hartwig.hmftools.esvee.assembly.alignment.Breakend;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;

import org.jfree.svg.SVGGraphics2D;

import htsjdk.samtools.CigarElement;
import j2html.tags.DomContent;

public class AssemblyVisualiser
{
    public record SegmentViewModel(String chromosome, BaseRegion viewRegion, BaseRegion refViewRegion, BaseSeqViewModel refViewModel, BaseSeqViewModel assemblyViewModel, boolean isInsert) {}

    private final AssemblyConfig mConfig;
    private final AssemblyAlignment mAssemblyAlignment;
    private final List<SegmentViewModel> mRefViewModel;

    public AssemblyVisualiser(final AssemblyConfig config, final AssemblyAlignment assemblyAlignment)
    {
        mConfig = config;
        mAssemblyAlignment = assemblyAlignment;
        mRefViewModel = getRefViewModel(config, assemblyAlignment);
    }

    public void writeVisFiles()
    {
        String info = mAssemblyAlignment.info();
        String filename = getFilename(info);

        CssBuilder readTableStyle = CssBuilder.EMPTY.borderSpacing(CssSize.ZERO);
        CssBuilder verticalSpacerStyle = CssBuilder.EMPTY.height(CssSize.em(1));

        DomContent verticalSpacer = div().withStyle(verticalSpacerStyle.toString());
        DomContent variantInfo = renderVariantInfo();

        List<DomContent> readTableRows = renderReadTable();
        DomContent readTable = div(styledTable(readTableRows, readTableStyle));
        String htmlStr = html(
                header(JQUERY_SCRIPT),
                body(
                        variantInfo,
                        verticalSpacer,
                        readTable,
                        getJavascript()).withStyle(BASE_FONT_STYLE.toString())).render();

        Path filePath = (new File(new File(mConfig.OutputDir, VIS_DIR), filename)).toPath();
        SV_LOGGER.debug("writing assembly vis file: {}", filePath.toString());

        try
        {
            Files.createDirectories(filePath.getParent());
            BufferedWriter outputWriter = createBufferedWriter(filePath.toString(), false);
            outputWriter.write(htmlStr);
            outputWriter.newLine();
            closeBufferedWriter(outputWriter);
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to write output file({}): {}", filePath.toString(), e.toString());
            System.exit(1);
        }
    }

    private DomContent renderVariantInfo()
    {
        CssBuilder strongStyle = CssBuilder.EMPTY.fontWeight("bold");

        StringJoiner infoJoiner = new StringJoiner(" ");
        infoJoiner.add(format("%s[normal: %s]", mConfig.TumorIds.get(0), mConfig.ReferenceIds.get(0)));
        infoJoiner.add(mAssemblyAlignment.toString());
        return div(infoJoiner.toString()).withStyle(strongStyle.toString());
    }

    private List<DomContent> renderReadTable()
    {
        CssBuilder lightGrayBgStyle = CssBuilder.EMPTY.backgroundColor(Color.LIGHT_GRAY);
        CssBuilder headerStyle = CssBuilder.EMPTY.backgroundColor(Color.LIGHT_GRAY).textAlign("right");

        List<DomContent> tableRows = Lists.newArrayList();

        // header row
        List<DomContent> headerCols = Lists.newArrayList();

        headerCols.add(td(""));
        headerCols.add(td(renderCoords()).withStyle(lightGrayBgStyle.toString()));
        DomContent headerRow = tr().with(headerCols);
        tableRows.add(headerRow);

        // chr label row
        DomContent chrLabelRow = tr(td("chr").withStyle(headerStyle.toString()), td(renderChomosomeLabels()));
        tableRows.add(chrLabelRow);

        // ref row
        DomContent refRow = tr(td("ref").withStyle(headerStyle.toString()), td(renderRef()));
        tableRows.add(refRow);

        // assembly row
        DomContent assemblyRow = tr(td("assembly").withStyle(headerStyle.toString()), td(renderAssembly()));
        tableRows.add(assemblyRow);

        // reads
        List<ReadViewModel> readViewModels = Lists.newArrayList();
        for(JunctionAssembly assembly : mAssemblyAlignment.assemblies())
        {
            for(SupportRead read : assembly.support())
            {
                if(read.isSupplementary())
                    continue;

                readViewModels.add(ReadViewModel.create(mRefViewModel, read, assembly));
            }
        }

        for(ReadViewModel readViewModel : readViewModels)
        {
            DomContent readEl = readViewModel.render();
            if(readEl == null)
                continue;

            tableRows.add(tr(td(""), td(readEl)));
        }

        return tableRows;
    }

    private record BreakendInfo(String chromosome, int pos, Orientation orient, List<CigarElement> cigarElements, String insertedBases, BaseRegion refRegion, String assemblySeq, String refSeq) {}

    private static BreakendInfo extractBreakendInfo(final RefGenomeInterface refGenome, final Breakend breakend, final String fullAssemblySeq)
    {
        AlignData alignment = breakend.segments().get(0).Alignment;
        List<CigarElement> cigarEls = alignment.cigarElements();
        String assemblySeq = fullAssemblySeq.substring(alignment.sequenceStart(), alignment.sequenceEnd() + 1);
        if(breakend.Orient == FORWARD)
        {
            if(cigarEls.get(cigarEls.size() - 1).getOperator().isClipping())
                cigarEls = cigarEls.subList(0, cigarEls.size() - 1);
        }
        else
        {
            if(cigarEls.get(0).getOperator().isClipping())
                cigarEls = cigarEls.subList(1, cigarEls.size());
        }

        int refLength = 0;
        for(CigarElement cigarEl : cigarEls)
        {
            if(cigarEl.getOperator().consumesReferenceBases())
                refLength += cigarEl.getLength();
        }

        int refStart;
        int refEnd;
        if(breakend.Orient == FORWARD)
        {
            refEnd = breakend.Position;
            refStart = refEnd - refLength + 1;
        }
        else
        {
            refStart = breakend.Position;
            refEnd = refStart + refLength - 1;
        }

        BaseRegion refRegion = new BaseRegion(refStart, refEnd);
        String refSeq = refGenome.getBaseString(breakend.Chromosome, refStart, refEnd);
        return new BreakendInfo(breakend.Chromosome, breakend.Position, breakend.Orient, cigarEls, breakend.InsertedBases, refRegion, assemblySeq, refSeq);
    }

    private static List<SegmentViewModel> getRefViewModel(final AssemblyConfig config, final AssemblyAlignment assemblyAlignment)
    {
        List<SegmentViewModel> refViewModel = Lists.newArrayList();

        String fullAssemblySeq = assemblyAlignment.fullSequence();
        int baseIdx = 0;
        for(int i = 0; i < assemblyAlignment.breakends().size(); i++)
        {
            BreakendInfo breakendInfo = extractBreakendInfo(config.RefGenome, assemblyAlignment.breakends().get(i), fullAssemblySeq);

            BaseSeqViewModel refSeqViewModel = BaseSeqViewModel.fromStr(breakendInfo.refSeq, baseIdx);
            BaseSeqViewModel assemblySeqViewModel = BaseSeqViewModel.fromStringWithCigar(breakendInfo.assemblySeq, breakendInfo.cigarElements, baseIdx);
            baseIdx += breakendInfo.refRegion.baseLength();

            BaseRegion viewRegion;
            BaseRegion refViewRegion;
            if(breakendInfo.orient == FORWARD)
            {
                int viewRegionEnd = refSeqViewModel.LastBasePos;
                int viewRegionStart = max(refSeqViewModel.FirstBasePos, viewRegionEnd - VIEW_REGION_SIZE + 1);
                viewRegion = new BaseRegion(viewRegionStart, viewRegionEnd);

                viewRegionEnd = breakendInfo.refRegion.end();
                viewRegionStart = max(breakendInfo.refRegion.start(), viewRegionEnd - VIEW_REGION_SIZE + 1);
                refViewRegion = new BaseRegion(viewRegionStart, viewRegionEnd);
            }
            else
            {
                int viewRegionStart = refSeqViewModel.FirstBasePos;
                int viewRegionEnd = min(refSeqViewModel.LastBasePos, viewRegionStart + VIEW_REGION_SIZE - 1);
                viewRegion = new BaseRegion(viewRegionStart, viewRegionEnd);

                viewRegionStart = breakendInfo.refRegion.start();
                viewRegionEnd = min(breakendInfo.refRegion.end(), viewRegionStart + VIEW_REGION_SIZE - 1);
                refViewRegion = new BaseRegion(viewRegionStart, viewRegionEnd);
            }

            refViewModel.add(new SegmentViewModel(breakendInfo.chromosome, viewRegion, refViewRegion, refSeqViewModel, assemblySeqViewModel, false));
            if(i < assemblyAlignment.breakends().size() - 1 && breakendInfo.insertedBases != null && !breakendInfo.insertedBases.isEmpty())
            {
                BaseSeqViewModel insertSeqViewModel = BaseSeqViewModel.fromStr(breakendInfo.insertedBases, baseIdx);
                baseIdx += breakendInfo.insertedBases.length();
                viewRegion = new BaseRegion(insertSeqViewModel.FirstBasePos, insertSeqViewModel.LastBasePos);
                refViewModel.add(new SegmentViewModel(null, viewRegion, null, null, insertSeqViewModel, true));
            }
        }

        return refViewModel;
    }

    private DomContent renderRef()
    {
        int totalBoxWidth = 0;
        for(SegmentViewModel refEl : mRefViewModel)
            totalBoxWidth += refEl.viewRegion.baseLength() + 2 * BOX_PADDING;

        SVGGraphics2D svgCanvas = new SVGGraphics2D(READ_HEIGHT_PX * totalBoxWidth, READ_HEIGHT_PX);
        AffineTransform initTransform = svgCanvas.getTransform();
        double xBoxOffset = 0.0d;
        for(SegmentViewModel refEl : mRefViewModel)
        {
            BaseSeqViewModel viewModel = refEl.refViewModel;
            BaseRegion viewRegion = refEl.viewRegion;
            if(viewModel == null)
            {
                xBoxOffset += viewRegion.baseLength() + 2 * BOX_PADDING;
                continue;
            }

            svgCanvas.setTransform(initTransform);
            renderBaseSeq(svgCanvas, new Point2D.Double(xBoxOffset, 0.0d), READ_HEIGHT_PX, viewRegion, viewModel, false, Maps.newHashMap(), null);
            xBoxOffset += viewRegion.baseLength() + 2 * BOX_PADDING;
        }

        return rawHtml(svgCanvas.getSVGElement());
    }

    private DomContent renderAssembly()
    {
        int totalBoxWidth = 0;
        for(SegmentViewModel refEl : mRefViewModel)
            totalBoxWidth += refEl.viewRegion.baseLength() + 2 * BOX_PADDING;

        SVGGraphics2D svgCanvas = new SVGGraphics2D(READ_HEIGHT_PX * totalBoxWidth, READ_HEIGHT_PX);
        AffineTransform initTransform = svgCanvas.getTransform();
        double xBoxOffset = 0.0d;
        for(SegmentViewModel refEl : mRefViewModel)
        {
            BaseSeqViewModel viewModel = refEl.assemblyViewModel;
            BaseRegion viewRegion = refEl.viewRegion;
            svgCanvas.setTransform(initTransform);
            renderBaseSeq(svgCanvas, new Point2D.Double(xBoxOffset, 0.0d), READ_HEIGHT_PX, viewRegion, viewModel, false, Maps.newHashMap(), refEl.refViewModel());
            xBoxOffset += viewRegion.baseLength() + 2 * BOX_PADDING;
        }

        return rawHtml(svgCanvas.getSVGElement());
    }

    private DomContent renderCoords()
    {
        double scalingFactor = READ_HEIGHT_PX / BASE_BOX_SIZE;
        double charWidth = getStringBounds(COORD_FONT, "9").getWidth();
        double maxStringWidth = 0.0d;
        int totalBoxWidth = 0;
        for(SegmentViewModel refEl : mRefViewModel)
        {
            if(refEl.isInsert)
            {
                totalBoxWidth += refEl.viewRegion.baseLength() + 2 * BOX_PADDING;
                continue;
            }

            BaseRegion refViewRegion = refEl.refViewRegion;
            totalBoxWidth += refViewRegion.baseLength() + 2 * BOX_PADDING;
            for(int i = refViewRegion.start(); i <= refViewRegion.end(); i++)
                maxStringWidth = max(maxStringWidth, getStringBounds(COORD_FONT, format(Locale.US, "%,d", i)).getWidth());
        }

        SVGGraphics2D svgCanvas = new SVGGraphics2D(
                scalingFactor * BASE_BOX_SIZE * totalBoxWidth,
                scalingFactor * (maxStringWidth + 2 * charWidth));
        AffineTransform initTransform = svgCanvas.getTransform();
        double xBoxOffset = 0.0d;
        Queue<JunctionAssembly> assemblyQueue = new ArrayDeque<>(mAssemblyAlignment.assemblies());
        for(SegmentViewModel refEl : mRefViewModel)
        {
            if(refEl.isInsert)
            {
                xBoxOffset += refEl.viewRegion.baseLength() + 2 * BOX_PADDING;
                continue;
            }

            BaseRegion refViewRegion = refEl.refViewRegion;
            Junction junction = assemblyQueue.poll().junction();
            int centerPosition = junction.Position;
            Point2D.Double canvasSize = new Point2D.Double(
                    scalingFactor * BASE_BOX_SIZE * (refViewRegion.baseLength() + 2 * BOX_PADDING), svgCanvas.getHeight());
            svgCanvas.setTransform(initTransform);
            SvgRender.renderCoords(svgCanvas, new Point2D.Double(xBoxOffset, 0.0d), canvasSize, READ_HEIGHT_PX, refViewRegion, centerPosition, DISPLAY_EVERY_NTH_COORD, maxStringWidth);
            xBoxOffset += refViewRegion.baseLength() + 2 * BOX_PADDING;
        }

        return rawHtml(svgCanvas.getSVGElement());
    }

    private DomContent renderChomosomeLabels()
    {
        List<ChrBaseRegion> regions = Lists.newArrayList();
        int xBoxOffset = 0;
        for(SegmentViewModel refEl : mRefViewModel)
        {
            BaseRegion viewRegion = refEl.viewRegion;
            int length = viewRegion.baseLength() + 2 * BOX_PADDING;
            if(refEl.isInsert)
            {
                xBoxOffset += length;
                continue;
            }

            regions.add(new ChrBaseRegion(refEl.chromosome, xBoxOffset, xBoxOffset + length - 1));
            xBoxOffset += length;
        }

        SVGGraphics2D svgCanvas = SvgRender.renderChrLabels(READ_HEIGHT_PX, regions);
        return rawHtml(svgCanvas.getSVGElement());
    }

    private String getFilename(final String info)
    {
        String sampleId = mConfig.TumorIds.get(0);
        return sampleId + "_" + info.replace(':', '_') + ".html";
    }
}
