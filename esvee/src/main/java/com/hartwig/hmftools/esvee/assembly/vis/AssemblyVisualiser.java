package com.hartwig.hmftools.esvee.assembly.vis;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.vis.BaseSeqViewModel.fromStr;
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
import java.util.List;
import java.util.Locale;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.vis.BaseSeqViewModel;
import com.hartwig.hmftools.common.vis.CssBuilder;
import com.hartwig.hmftools.common.vis.CssSize;
import com.hartwig.hmftools.common.vis.SvgRender;
import com.hartwig.hmftools.esvee.assembly.AssemblyConfig;
import com.hartwig.hmftools.esvee.assembly.alignment.AlignData;
import com.hartwig.hmftools.esvee.assembly.alignment.AssemblyAlignment;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;

import org.jfree.svg.SVGGraphics2D;

import j2html.tags.DomContent;

public class AssemblyVisualiser
{
    public record RefSegmentViewModel(JunctionAssembly assembly, BaseRegion viewRegion, BaseSeqViewModel viewModel) {}

    private final AssemblyConfig mConfig;
    private final AssemblyAlignment mAssemblyAlignment;

    private final List<RefSegmentViewModel> mRefViewModel;

    public AssemblyVisualiser(final AssemblyConfig config, final AssemblyAlignment assemblyAlignment)
    {
        mConfig = config;
        mAssemblyAlignment = assemblyAlignment;
        mRefViewModel = getRefViewModel(assemblyAlignment);
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

        // reads
        List<ReadViewModel> readViewModels = Lists.newArrayList();
        for(JunctionAssembly assembly : mAssemblyAlignment.assemblies())
        {
            for(SupportRead read : assembly.support())
            {
                if(read.isSupplementary())
                    continue;

                readViewModels.add(ReadViewModel.create(mRefViewModel, assembly, read));
            }
        }

        for(ReadViewModel readViewModel : readViewModels)
        {
            List<DomContent> cols = Lists.newArrayList();
            cols.add(td(""));
            cols.add(td(readViewModel.render()));
            tableRows.add(tr().with(cols));
        }

        return tableRows;
    }

    private static List<RefSegmentViewModel> getRefViewModel(final AssemblyAlignment assemblyAlignment)
    {
        String refSeq = assemblyAlignment.fullSequence();
        List<AlignData> segments = assemblyAlignment.breakends().stream()
                .flatMap(x -> x.segments().stream())
                .map(x -> x.Alignment)
                .toList();

        List<RefSegmentViewModel> refViewModel = Lists.newArrayList();
        List<Integer> linkIndices = assemblyAlignment.linkIndices();
        int refIdx = 0;
        for(int i = 0; i < segments.size(); i++)
        {
            AlignData segment = segments.get(i);

            String chromosome = segment.chromosome();
            int posStart = segment.positionStart();
            int posEnd = segment.positionEnd();
            int length = posEnd - posStart + 1;

            // find associated JunctionAssembly
            JunctionAssembly junctionAssembly = null;
            for(JunctionAssembly assembly : assemblyAlignment.assemblies())
            {
                Junction junction = assembly.junction();
                if(!junction.Chromosome.equals(chromosome))
                    continue;

                if(junction.Orient == FORWARD && posEnd != junction.Position)
                    continue;

                if(junction.Orient == REVERSE && posStart != junction.Position)
                    continue;

                junctionAssembly = assembly;
                break;
            }

            BaseSeqViewModel segmentViewModel = fromStr(refSeq.substring(refIdx, refIdx + length), chromosome, posStart);
            BaseRegion viewRegion;
            if(junctionAssembly.junction().Orient == FORWARD)
            {
                int viewRegionEnd = segmentViewModel.LastBasePos;
                int viewRegionStart = max(segmentViewModel.FirstBasePos, viewRegionEnd - VIEW_REGION_SIZE + 1);
                viewRegion = new BaseRegion(viewRegionStart, viewRegionEnd);
            }
            else
            {
                int viewRegionStart = segmentViewModel.FirstBasePos;
                int viewRegionEnd = min(segmentViewModel.LastBasePos, viewRegionStart + VIEW_REGION_SIZE - 1);
                viewRegion = new BaseRegion(viewRegionStart, viewRegionEnd);
            }

            refViewModel.add(new RefSegmentViewModel(junctionAssembly, viewRegion, segmentViewModel));
            if(i < segments.size() - 1)
                refIdx = linkIndices.get(i);
        }

        return refViewModel;
    }

    private DomContent renderRef()
    {
        int totalBoxWidth = 0;
        for(RefSegmentViewModel segmentViewModel : mRefViewModel)
            totalBoxWidth += segmentViewModel.viewRegion.baseLength() + 2 * BOX_PADDING;

        SVGGraphics2D svgCanvas = new SVGGraphics2D(READ_HEIGHT_PX * totalBoxWidth, READ_HEIGHT_PX);
        AffineTransform initTransform = svgCanvas.getTransform();
        double xBoxOffset = 0.0d;
        for(RefSegmentViewModel segmentViewModel : mRefViewModel)
        {
            BaseSeqViewModel viewModel = segmentViewModel.viewModel;
            BaseRegion viewRegion = segmentViewModel.viewRegion;
            svgCanvas.setTransform(initTransform);
            renderBaseSeq(svgCanvas, new Point2D.Double(xBoxOffset, 0.0d), READ_HEIGHT_PX, viewRegion, viewModel, false, Maps.newHashMap(), null);
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
        for(RefSegmentViewModel segmentViewModel : mRefViewModel)
        {
            BaseRegion viewRegion = segmentViewModel.viewRegion;
            totalBoxWidth += viewRegion.baseLength() + 2 * BOX_PADDING;
            for(int i = viewRegion.start(); i <= viewRegion.end(); i++)
                maxStringWidth = max(maxStringWidth, getStringBounds(COORD_FONT, format(Locale.US, "%,d", i)).getWidth());
        }

        SVGGraphics2D svgCanvas = new SVGGraphics2D(
                scalingFactor * BASE_BOX_SIZE * totalBoxWidth,
                scalingFactor * (maxStringWidth + 2 * charWidth));
        AffineTransform initTransform = svgCanvas.getTransform();
        double xBoxOffset = 0.0d;
        for(int i = 0; i < mRefViewModel.size(); i++)
        {
            RefSegmentViewModel segmentViewModel = mRefViewModel.get(i);
            Junction junction = mAssemblyAlignment.assemblies().get(i).junction();

            BaseRegion viewRegion = segmentViewModel.viewRegion;
            svgCanvas.setTransform(initTransform);

            int centerPosition = junction.Position;
            Point2D.Double canvasSize = new Point2D.Double(
                    scalingFactor * BASE_BOX_SIZE * (viewRegion.baseLength() + 2 * BOX_PADDING), svgCanvas.getHeight());
            SvgRender.renderCoords(svgCanvas, new Point2D.Double(xBoxOffset, 0.0d), canvasSize, READ_HEIGHT_PX, viewRegion, centerPosition, DISPLAY_EVERY_NTH_COORD, maxStringWidth);
            xBoxOffset += viewRegion.baseLength() + 2 * BOX_PADDING;
        }

        return rawHtml(svgCanvas.getSVGElement());
    }

    private DomContent renderChomosomeLabels()
    {
        List<ChrBaseRegion> regions = Lists.newArrayList();
        int xBoxOffset = 0;
        for(RefSegmentViewModel segmentViewModel : mRefViewModel)
        {
            BaseSeqViewModel viewModel = segmentViewModel.viewModel;
            BaseRegion viewRegion = segmentViewModel.viewRegion;
            int length = viewRegion.baseLength() + 2 * BOX_PADDING;
            regions.add(new ChrBaseRegion(viewModel.Chromosome, xBoxOffset, xBoxOffset + length - 1));
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
