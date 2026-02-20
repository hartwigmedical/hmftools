package com.hartwig.hmftools.esvee.assembly.vis;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.codon.Nucleotides.reverseComplementBases;
import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.sv.SvVcfTags.ASM_ID;
import static com.hartwig.hmftools.common.sv.SvVcfTags.ASM_INFO;
import static com.hartwig.hmftools.common.sv.SvVcfTags.ASM_LENGTH;
import static com.hartwig.hmftools.common.sv.SvVcfTags.AVG_FRAG_LENGTH;
import static com.hartwig.hmftools.common.sv.SvVcfTags.BE_ASM_ORIENT;
import static com.hartwig.hmftools.common.sv.SvVcfTags.BE_ASM_POS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.BE_ORIENT;
import static com.hartwig.hmftools.common.sv.SvVcfTags.DISC_FRAGS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.MATE_ID;
import static com.hartwig.hmftools.common.sv.SvVcfTags.MAX_LOCAL_REPEAT;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SEG_ALIGN_LENGTH;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SEG_MAPQ;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SEG_REPEAT_LENGTH;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SEG_SCORE;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SPLIT_FRAGS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SV_ID;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SV_TYPE;
import static com.hartwig.hmftools.common.sv.SvVcfTags.THREE_PRIME_RANGE;
import static com.hartwig.hmftools.common.sv.SvVcfTags.TOTAL_FRAGS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.UNIQUE_FRAG_POSITIONS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.VCF_ITEM_DELIM;
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
import static com.hartwig.hmftools.common.vis.SvgUtil.Alignment.LEFT;
import static com.hartwig.hmftools.common.vis.SvgUtil.Alignment.RIGHT;
import static com.hartwig.hmftools.common.vis.SvgUtil.getStringBounds;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.BREAKEND_REQ_VALID_FRAGMENT_LENGTH_PERC;
import static com.hartwig.hmftools.esvee.assembly.output.VcfWriter.buildGenotype;
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
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.NavigableMap;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.vis.BaseSeqViewModel;
import com.hartwig.hmftools.common.vis.CssBuilder;
import com.hartwig.hmftools.common.vis.CssSize;
import com.hartwig.hmftools.common.vis.SvgRender;
import com.hartwig.hmftools.common.vis.SvgRender.ChrLabel;
import com.hartwig.hmftools.common.vis.SvgUtil.Alignment;
import com.hartwig.hmftools.esvee.assembly.AssemblyConfig;
import com.hartwig.hmftools.esvee.assembly.alignment.AlignData;
import com.hartwig.hmftools.esvee.assembly.alignment.AssemblyAlignment;
import com.hartwig.hmftools.esvee.assembly.alignment.Breakend;
import com.hartwig.hmftools.esvee.assembly.alignment.BreakendSegment;
import com.hartwig.hmftools.esvee.assembly.alignment.BreakendSupport;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;

import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.Nullable;
import org.jfree.svg.SVGGraphics2D;

import htsjdk.samtools.CigarElement;
import htsjdk.variant.variantcontext.Genotype;
import j2html.tags.DomContent;

public class AssemblyVisualiser
{
    private static final String MISSING_FIELD = "missing";

    public record SegmentViewModel(Integer breakendId, String chromosome, @Nullable Integer position, BaseRegion refRegion,
                                   BaseRegion viewRegion,
                                   BaseRegion refViewRegion, BaseSeqViewModel refViewModel, boolean isRefReversed,
                                   boolean isAssemblyReversed,
                                   BaseSeqViewModel assemblyViewModel, boolean isInsert, int leftDelLength, BaseRegion assemblyRegion) {}

    public record PairedSegmentViewModel(List<SegmentViewModel> viewModels, int readStartOffset, int prevBoxWidth) {}

    private final AssemblyConfig mConfig;
    private final AssemblyAlignment mAssemblyAlignment;
    private final List<PairedSegmentViewModel> mRefViewModel;

    private final String mFullAssemblyStr;

    private final List<String> mSampleNames;
    private final Map<String, Integer> mSampleNameIndex;
    private final Set<String> mTumorIds;

    public AssemblyVisualiser(final AssemblyConfig config, final AssemblyAlignment assemblyAlignment)
    {
        mConfig = config;
        mAssemblyAlignment = assemblyAlignment;
        mRefViewModel = getRefViewModel(config, assemblyAlignment);

        mFullAssemblyStr = assemblyAlignment.fullSequence();

        mTumorIds = Sets.newHashSet(mConfig.TumorIds);
        mSampleNames = Lists.newArrayList(mConfig.ReferenceIds);
        mSampleNames.addAll(mConfig.TumorIds);

        mSampleNameIndex = Maps.newHashMap();
        List<String> configSampleIds = mConfig.combinedSampleIds();
        for(int i = 0; i < configSampleIds.size(); ++i)
        {
            String sampleId = configSampleIds.get(i);
            mSampleNameIndex.put(sampleId, i);
        }
    }

    public void writeVisFiles()
    {
        String info = mAssemblyAlignment.info();
        String filename = getFilename(info);

        CssBuilder readTableStyle = CssBuilder.EMPTY.borderSpacing(CssSize.ZERO);
        CssBuilder verticalSpacerStyle = CssBuilder.EMPTY.height(CssSize.em(1));

        DomContent verticalSpacer = div().withStyle(verticalSpacerStyle.toString());
        DomContent variantInfo = renderVariantInfo();

        // TODO(mkcmkc): remove try/catch
        List<DomContent> readTableRows;
        try
        {
            readTableRows = renderReadTable();
        }
        catch(Exception e)
        {
            e.printStackTrace();
            return;
        }

        List<DomContent> bodyElements = Lists.newArrayList();
        bodyElements.add(variantInfo);
        bodyElements.add(verticalSpacer);

        DomContent sampleInfo = styledTable(List.of(tr(
                td(renderVariantAndSampleSupportTable(mAssemblyAlignment.breakends()))
                        .withStyle(CssBuilder.EMPTY.verticalAlign("top").toString()),
                td(renderPropertiesTable(mAssemblyAlignment.breakends()))
                        .withStyle(CssBuilder.EMPTY.verticalAlign("top").toString()),
                td(renderAssemblyInfo())
                        .withStyle(CssBuilder.EMPTY.verticalAlign("top").toString())
        )), CssBuilder.EMPTY);

        bodyElements.add(sampleInfo);
        bodyElements.add(verticalSpacer);

        DomContent readTable = div(styledTable(readTableRows, readTableStyle));
        bodyElements.add(readTable);
        bodyElements.add(getJavascript());

        String htmlStr = html(
                header(JQUERY_SCRIPT),
                body().with(bodyElements).withStyle(BASE_FONT_STYLE.toString())).render();

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
        return div(format("%s[normal: %s]", mConfig.TumorIds.get(0), mConfig.ReferenceIds.get(0))).withStyle(strongStyle.toString());
    }

    private static void renderInfoTableSection(final List<DomContent> rows, final String sectionHeader,
            final LinkedHashMap<String, List<String>> info, final CssBuilder leftHeaderStyle, final CssBuilder firstRowCellStyle,
            final CssBuilder cellStyle)
    {
        boolean firstRow = true;
        for(Map.Entry<String, List<String>> entry : info.entrySet())
        {
            List<DomContent> rowElements = Lists.newArrayList();
            if(firstRow)
                rowElements.add(td(sectionHeader).attr("rowspan", info.size()).withStyle(leftHeaderStyle.toString()));

            CssBuilder elementStyle;
            if(firstRow)
                elementStyle = firstRowCellStyle.textAlign("right");
            else
                elementStyle = cellStyle.textAlign("right");

            rowElements.add(td(entry.getKey()).withStyle(elementStyle.toString()));
            for(String value : entry.getValue())
                rowElements.add(td(value).withStyle(elementStyle.toString()));

            rows.add(tr().with(rowElements));
            firstRow = false;
        }
    }

    private DomContent renderVariantAndSampleSupportTable(final List<Breakend> breakends)
    {
        CssBuilder borderStyle = CssBuilder.EMPTY.border(CssSize.px(1.0), "solid", Color.BLACK).borderCollapse("collapse");
        CssBuilder firstRowBorderStyle = CssBuilder.EMPTY.borderLeft(CssSize.px(1.0), "solid", Color.BLACK)
                .borderRight(CssSize.px(1.0), "solid", Color.BLACK)
                .borderTop(CssSize.px(2.0), "solid", Color.BLACK)
                .borderBottom(CssSize.px(1.0), "solid", Color.BLACK)
                .borderCollapse("collapse");
        CssBuilder tableStyle = borderStyle.marginRight(CssSize.px(10.0));
        CssBuilder headerStyle = CssBuilder.EMPTY.fontWeight("bold").textAlign("center").backgroundColor(Color.LIGHT_GRAY);
        CssBuilder cellStyle = CssBuilder.EMPTY.paddingLeft(CssSize.em(0.5)).paddingRight(CssSize.em(0.5)).merge(borderStyle);
        CssBuilder firstRowCellStyle = CssBuilder.EMPTY.paddingLeft(CssSize.em(0.5))
                .paddingRight(CssSize.em(0.5))
                .merge(firstRowBorderStyle);
        CssBuilder leftHeaderStyle = CssBuilder.EMPTY.backgroundColor(Color.LIGHT_GRAY)
                .writingMode("vertical-rl")
                .transform("rotate(180deg)")
                .fontWeight("bold")
                .textAlign("center")
                .merge(firstRowBorderStyle);

        LinkedHashMap<String, List<String>> variantInfo = new LinkedHashMap<>();
        variantInfo.put("Chromosome", breakends.stream().map(x -> x.Chromosome).toList());
        variantInfo.put("Position", breakends.stream().map(x -> String.valueOf(x.Position)).toList());
        variantInfo.put("Qual", breakends.stream().map(x -> String.valueOf(x.calcSvQual())).toList());
        variantInfo.put(SV_TYPE, breakends.stream().map(x -> x.svType().name()).toList());

        LinkedHashMap<String, List<String>> sampleSupportInfo = new LinkedHashMap<>();
        List<List<Genotype>> genotypes = Lists.newArrayList();
        for(Breakend breakend : breakends)
        {
            List<Genotype> breakendGenotypes = Lists.newArrayList();
            for(String sampleId : mSampleNames)
            {
                int sampleSupportIndex = mSampleNameIndex.get(sampleId);
                BreakendSupport breakendSupport = breakend.sampleSupport().get(sampleSupportIndex);
                breakendGenotypes.add(buildGenotype(sampleId, breakendSupport));
            }

            genotypes.add(breakendGenotypes);
        }

        List<Integer> breakendTotalDiscordantFrags = Lists.newArrayList();
        List<Integer> breakendTotalSplitFrags = Lists.newArrayList();
        for(Breakend breakend : breakends)
        {
            int totalSplitFrags = 0;
            int totalDiscFrags = 0;
            for(String sampleId : mSampleNames)
            {
                int sampleSupportIndex = mSampleNameIndex.get(sampleId);
                BreakendSupport breakendSupport = breakend.sampleSupport().get(sampleSupportIndex);
                totalSplitFrags += breakendSupport.SplitFragments;
                totalDiscFrags += breakendSupport.DiscordantFragments;
            }

            breakendTotalDiscordantFrags.add(totalDiscFrags);
            breakendTotalSplitFrags.add(totalSplitFrags);
        }

        for(int i = 0; i < breakends.size(); i++)
        {
            List<Genotype> breakendGenotypes = genotypes.get(i);
            int totalDiscordantFrags = breakendTotalDiscordantFrags.get(i);
            int totalSplitFrags = breakendTotalSplitFrags.get(i);
            Genotype tumorGenotype = null;
            for(Genotype genotype : breakendGenotypes)
            {
                if(!mTumorIds.contains(genotype.getSampleName()))
                    continue;

                tumorGenotype = genotype;
                break;
            }

            if(tumorGenotype == null)
            {
                sampleSupportInfo.computeIfAbsent("AD", k -> Lists.newArrayList()).add(MISSING_FIELD);
                sampleSupportInfo.computeIfAbsent(TOTAL_FRAGS, k -> Lists.newArrayList()).add(MISSING_FIELD);
                sampleSupportInfo.computeIfAbsent(DISC_FRAGS, k -> Lists.newArrayList()).add(MISSING_FIELD);
                sampleSupportInfo.computeIfAbsent("DP", k -> Lists.newArrayList()).add(MISSING_FIELD);
                sampleSupportInfo.computeIfAbsent("REF", k -> Lists.newArrayList()).add(MISSING_FIELD);
                sampleSupportInfo.computeIfAbsent("REFPAIR", k -> Lists.newArrayList()).add(MISSING_FIELD);
                sampleSupportInfo.computeIfAbsent("SB", k -> Lists.newArrayList()).add(MISSING_FIELD);
                sampleSupportInfo.computeIfAbsent(SPLIT_FRAGS, k -> Lists.newArrayList()).add(MISSING_FIELD);
            }
            else
            {
                String adString = Arrays.stream(tumorGenotype.getAD()).mapToObj(String::valueOf).collect(Collectors.joining(","));
                sampleSupportInfo.computeIfAbsent("AD", k -> Lists.newArrayList()).add(adString);
                sampleSupportInfo.computeIfAbsent(TOTAL_FRAGS, k -> Lists.newArrayList())
                        .add(String.valueOf(totalDiscordantFrags + totalSplitFrags));
                sampleSupportInfo.computeIfAbsent(DISC_FRAGS, k -> Lists.newArrayList()).add(String.valueOf(totalDiscordantFrags));
                sampleSupportInfo.computeIfAbsent("DP", k -> Lists.newArrayList()).add(String.valueOf(tumorGenotype.getDP()));

                // TODO(mkcmkc): These fields are not accessible from Assembly application.
                sampleSupportInfo.computeIfAbsent("REF", k -> Lists.newArrayList()).add("TODO");
                sampleSupportInfo.computeIfAbsent("REFPAIR", k -> Lists.newArrayList()).add("TODO");

                sampleSupportInfo.computeIfAbsent("SB", k -> Lists.newArrayList())
                        .add(format("%.3f", tumorGenotype.getExtendedAttribute("SB")));
                sampleSupportInfo.computeIfAbsent(SPLIT_FRAGS, k -> Lists.newArrayList()).add(String.valueOf(totalSplitFrags));
            }
        }

        List<DomContent> rows = Lists.newArrayList();

        List<DomContent> headerElements = Lists.newArrayList();
        headerElements.add(td().withStyle(cellStyle.merge(headerStyle).toString()));
        headerElements.add(td().withStyle(cellStyle.merge(headerStyle).toString()));
        for(int i = 0; i < breakends.size(); i++)
            headerElements.add(td(String.valueOf(i)).withStyle(cellStyle.merge(headerStyle).toString()));

        rows.add(tr().with(headerElements));

        renderInfoTableSection(rows, "Variant", variantInfo, leftHeaderStyle, firstRowCellStyle, cellStyle);
        renderInfoTableSection(rows, "Sample support", sampleSupportInfo, leftHeaderStyle, firstRowCellStyle, cellStyle);

        DomContent table = styledTable(rows, tableStyle);
        return div(table);
    }

    private static DomContent renderPropertiesTable(final List<Breakend> breakends)
    {
        CssBuilder borderStyle = CssBuilder.EMPTY.border(CssSize.px(1.0), "solid", Color.BLACK).borderCollapse("collapse");
        CssBuilder firstRowBorderStyle = CssBuilder.EMPTY.borderLeft(CssSize.px(1.0), "solid", Color.BLACK)
                .borderRight(CssSize.px(1.0), "solid", Color.BLACK)
                .borderTop(CssSize.px(2.0), "solid", Color.BLACK)
                .borderBottom(CssSize.px(1.0), "solid", Color.BLACK)
                .borderCollapse("collapse");
        CssBuilder tableStyle = borderStyle.marginRight(CssSize.px(10.0));
        CssBuilder headerStyle = CssBuilder.EMPTY.fontWeight("bold").textAlign("center").backgroundColor(Color.LIGHT_GRAY);
        CssBuilder cellStyle = CssBuilder.EMPTY.paddingLeft(CssSize.em(0.5)).paddingRight(CssSize.em(0.5)).merge(borderStyle);
        CssBuilder firstRowCellStyle = CssBuilder.EMPTY.paddingLeft(CssSize.em(0.5))
                .paddingRight(CssSize.em(0.5))
                .merge(firstRowBorderStyle);
        CssBuilder leftHeaderStyle = CssBuilder.EMPTY.backgroundColor(Color.LIGHT_GRAY)
                .writingMode("vertical-rl")
                .transform("rotate(180deg)")
                .fontWeight("bold")
                .textAlign("center")
                .merge(firstRowBorderStyle);

        LinkedHashMap<String, List<String>> idInfo = new LinkedHashMap<>();
        idInfo.put("Id", breakends.stream().map(x -> String.valueOf(x.id())).toList());

        for(Breakend breakend : breakends)
        {
            int otherBreakendId = breakend.otherBreakend().id();
            idInfo.computeIfAbsent(MATE_ID, k -> Lists.newArrayList()).add(String.valueOf(otherBreakendId));

            // always write lower first
            int lowerBreakendId = min(otherBreakendId, breakend.id());
            int upperBreakendId = max(otherBreakendId, breakend.id());
            String svId = format("%d_%d", lowerBreakendId, upperBreakendId);
            idInfo.computeIfAbsent(SV_ID, k -> Lists.newArrayList()).add(svId);
        }

        LinkedHashMap<String, List<String>> segmentPropertiesInfo = new LinkedHashMap<>();
        for(Breakend breakend : breakends)
        {
            List<BreakendSegment> segments = breakend.segments();
            segmentPropertiesInfo.computeIfAbsent(BE_ASM_POS, k -> Lists.newArrayList()).add(
                    segments.stream().map(x -> String.valueOf(x.SequenceIndex)).collect(Collectors.joining(VCF_ITEM_DELIM)));

            segmentPropertiesInfo.computeIfAbsent(BE_ORIENT, k -> Lists.newArrayList()).add(breakend.Orient.name());
            segmentPropertiesInfo.computeIfAbsent(BE_ASM_ORIENT, k -> Lists.newArrayList()).add(
                    segments.stream().map(x -> x.Alignment.orientation().name()).collect(Collectors.joining(VCF_ITEM_DELIM)));

            segmentPropertiesInfo.computeIfAbsent(SEG_ALIGN_LENGTH, k -> Lists.newArrayList())
                    .add(segments.stream()
                            .map(x -> String.valueOf(x.Alignment.alignedBases()))
                            .collect(Collectors.joining(VCF_ITEM_DELIM)));
            segmentPropertiesInfo.computeIfAbsent(SEG_MAPQ, k -> Lists.newArrayList())
                    .add(String.valueOf(segments.stream().mapToInt(x -> x.Alignment.mapQual()).max().orElse(0)));
            segmentPropertiesInfo.computeIfAbsent(SEG_SCORE, k -> Lists.newArrayList())
                    .add(String.valueOf(segments.stream().mapToInt(x -> x.Alignment.score()).max().orElse(0)));
            segmentPropertiesInfo.computeIfAbsent(SEG_REPEAT_LENGTH, k -> Lists.newArrayList())
                    .add(String.valueOf(segments.stream().mapToInt(x -> x.Alignment.adjustedAlignment()).max().orElse(0)));
        }

        LinkedHashMap<String, List<String>> siteQualInfo = new LinkedHashMap<>();
        for(Breakend breakend : breakends)
        {
            if(breakend.isSingle() && breakend.validFragmentLengthPercent() < BREAKEND_REQ_VALID_FRAGMENT_LENGTH_PERC)
                siteQualInfo.computeIfAbsent(AVG_FRAG_LENGTH, k -> Lists.newArrayList()).add("0");
            else
                siteQualInfo.computeIfAbsent(AVG_FRAG_LENGTH, k -> Lists.newArrayList())
                        .add(String.valueOf(breakend.averageFragmentLength()));

            siteQualInfo.computeIfAbsent(UNIQUE_FRAG_POSITIONS, k -> Lists.newArrayList())
                    .add(String.valueOf(breakend.uniqueFragmentPositionCount()));
            siteQualInfo.computeIfAbsent(THREE_PRIME_RANGE, k -> Lists.newArrayList())
                    .add(Arrays.stream(breakend.readOrientationRange()).mapToObj(String::valueOf).collect(Collectors.joining(",")));
            siteQualInfo.computeIfAbsent(MAX_LOCAL_REPEAT, k -> Lists.newArrayList()).add(String.valueOf(breakend.maxLocalRepeat()));
        }

        List<DomContent> rows = Lists.newArrayList();

        List<DomContent> headerElements = Lists.newArrayList();
        headerElements.add(td().withStyle(cellStyle.merge(headerStyle).toString()));
        headerElements.add(td().withStyle(cellStyle.merge(headerStyle).toString()));
        for(int i = 0; i < breakends.size(); i++)
            headerElements.add(td(String.valueOf(i)).withStyle(cellStyle.merge(headerStyle).toString()));

        rows.add(tr().with(headerElements));

        renderInfoTableSection(rows, "Id", idInfo, leftHeaderStyle, firstRowCellStyle, cellStyle);
        renderInfoTableSection(rows, "Segment properties", segmentPropertiesInfo, leftHeaderStyle, firstRowCellStyle, cellStyle);
        renderInfoTableSection(rows, "Site Qual", siteQualInfo, leftHeaderStyle, firstRowCellStyle, cellStyle);

        DomContent table = styledTable(rows, tableStyle);
        return div(table);
    }

    private DomContent renderAssemblyInfo()
    {
        CssBuilder borderStyle = CssBuilder.EMPTY.border(CssSize.px(1.0), "solid", Color.BLACK).borderCollapse("collapse");
        CssBuilder tableStyle = borderStyle.marginRight(CssSize.px(10.0));
        CssBuilder cellStyle = CssBuilder.EMPTY.paddingLeft(CssSize.em(0.5)).paddingRight(CssSize.em(0.5)).merge(borderStyle);

        LinkedHashMap<String, String> info = new LinkedHashMap<>();
        info.put(ASM_INFO, mAssemblyAlignment.info(10));
        info.put(ASM_ID, String.valueOf(mAssemblyAlignment.id()));
        info.put(ASM_LENGTH, String.valueOf(mAssemblyAlignment.fullSequenceLength()));

        List<DomContent> infoRows = Lists.newArrayList();
        for(Map.Entry<String, String> entry : info.entrySet())
            infoRows.add(tr(
                    td(entry.getKey()).withStyle(cellStyle.textAlign("right").toString()),
                    td(entry.getValue()).withStyle(cellStyle.textAlign("right").toString())));

        return div(styledTable(infoRows, tableStyle));
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

                // TODO(mkcmkc): remove try
                try
                {
                    for(PairedSegmentViewModel pairedRefEl : mRefViewModel)
                    {
                        ReadViewModel readViewModel = ReadViewModel.create(pairedRefEl, read, assembly, mFullAssemblyStr);
                        if(readViewModel == null)
                            continue;

                        readViewModels.add(readViewModel);
                    }
                }
                catch(Exception e)
                {
                    e.printStackTrace();
                    throw e;
                }
            }
        }

        Collections.sort(readViewModels, Comparator.comparingInt(ReadViewModel::firstBaseIndex));
        for(ReadViewModel readViewModel : readViewModels)
        {
            DomContent readEl = readViewModel.render();
            if(readEl == null)
                continue;

            tableRows.add(tr(td(""), td(readEl)));
        }

        return tableRows;
    }

    private record BreakendInfo(Integer id, String chromosome, List<CigarElement> cigarElements, String insertedBases,
                                BaseRegion refRegion, String assemblySeq, String refSeq, boolean isRefReversed, boolean isAssemblyReversed,
                                BaseRegion assemblyRegion,
                                int leftDelLength) {}

    private static BreakendInfo extractBreakendInfo(final RefGenomeInterface refGenome, final PairedBreakendInfo pairedBreakend,
            final String fullAssemblySeq, @Nullable final Integer lastSeqEnd)
    {
        Breakend breakend = pairedBreakend.Breakend;
        AlignData alignment = breakend.segments().get(0).Alignment;
        List<CigarElement> cigarEls = alignment.cigarElements();
        int leftDelLength = 0;
        if(lastSeqEnd != null && lastSeqEnd >= alignment.sequenceStart())
            leftDelLength = lastSeqEnd - alignment.sequenceStart() + 1;

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

        boolean reverseAssemblySeq = pairedBreakend.IsAssemblyReversed;
        if(pairedBreakend.IsRefReversed)
        {
            Collections.reverse(cigarEls);
            reverseAssemblySeq = !reverseAssemblySeq;
        }

        if(reverseAssemblySeq)
        {
            assemblySeq = reverseComplementBases(assemblySeq);
            Collections.reverse(cigarEls);
        }

        String chromosome = alignment.refLocation().chromosome();
        BaseRegion refRegion = alignment.refLocation().baseRegion();
        String refSeq = refGenome.getBaseString(chromosome, refRegion.start(), refRegion.end());
        if(pairedBreakend.IsRefReversed)
            refSeq = reverseComplementBases(refSeq);

        BaseRegion assemblyRegion = new BaseRegion(alignment.sequenceStart(), alignment.sequenceEnd());
        return new BreakendInfo(breakend.id(), chromosome, cigarEls, breakend.InsertedBases, refRegion, assemblySeq, refSeq, pairedBreakend.IsRefReversed, reverseAssemblySeq, assemblyRegion, leftDelLength);
    }

    private record PairedRefViewModel(List<SegmentViewModel> refViewModelPair, int baseIdx, int boxWidth) {}

    private static PairedRefViewModel getPairedRefViewModel(final AssemblyConfig config, final int baseIdxStart, final String fullAssemblySeq,
            final Pair<PairedBreakendInfo, PairedBreakendInfo> breakendPair)
    {
        List<SegmentViewModel> refViewModelPair = Lists.newArrayList();
        Integer lastSequenceEnd = null;
        int baseIdx = baseIdxStart;
        int boxWidth = 0;
        for(int i = 0; i < 2; i++)
        {
            PairedBreakendInfo pairedBreakend = i == 0 ? breakendPair.getLeft() : breakendPair.getRight();
            BreakendInfo breakendInfo = extractBreakendInfo(config.RefGenome, pairedBreakend, fullAssemblySeq, lastSequenceEnd);
            lastSequenceEnd = breakendInfo.assemblyRegion.end();

            BaseSeqViewModel refSeqViewModel = BaseSeqViewModel.fromStr(breakendInfo.refSeq, baseIdx);
            BaseSeqViewModel assemblySeqViewModel = BaseSeqViewModel.fromStringWithCigar(
                    breakendInfo.assemblySeq, breakendInfo.cigarElements, baseIdx);

            BaseRegion viewRegion;
            if(i == 0)
            {
                int viewRegionEnd = baseIdx + breakendInfo.refRegion.baseLength() - 1;
                int viewRegionStart = max(baseIdx, viewRegionEnd - VIEW_REGION_SIZE + 1);
                viewRegion = new BaseRegion(viewRegionStart, viewRegionEnd);
            }
            else
            {
                int viewRegionStart = baseIdx;
                int viewRegionEnd = min(baseIdx + breakendInfo.refRegion.baseLength() - 1, viewRegionStart + VIEW_REGION_SIZE - 1);
                viewRegion = new BaseRegion(viewRegionStart, viewRegionEnd);
            }

            baseIdx += breakendInfo.refRegion.baseLength();

            boolean isViewRegionReversed = breakendInfo.isRefReversed;
            if(i != 0)
                isViewRegionReversed = !isViewRegionReversed;

            BaseRegion refViewRegion;
            int pos;
            if(!isViewRegionReversed)
            {
                int viewRegionEnd = breakendInfo.refRegion.end();
                int viewRegionStart = max(breakendInfo.refRegion.start(), viewRegionEnd - VIEW_REGION_SIZE + 1);
                refViewRegion = new BaseRegion(viewRegionStart, viewRegionEnd);
                pos = viewRegionEnd;
            }
            else
            {
                int viewRegionStart = breakendInfo.refRegion.start();
                int viewRegionEnd = min(breakendInfo.refRegion.end(), viewRegionStart + VIEW_REGION_SIZE - 1);
                refViewRegion = new BaseRegion(viewRegionStart, viewRegionEnd);
                pos = viewRegionStart;
            }

            refViewModelPair.add(new SegmentViewModel(breakendInfo.id, breakendInfo.chromosome, pos, breakendInfo.refRegion,
                    viewRegion, refViewRegion, refSeqViewModel, pairedBreakend.IsRefReversed, breakendInfo.isAssemblyReversed, assemblySeqViewModel, false, breakendInfo.leftDelLength, breakendInfo.assemblyRegion));
            boxWidth += refViewRegion.baseLength() + 2 * BOX_PADDING;
            if(i == 0 && breakendInfo.insertedBases != null && !breakendInfo.insertedBases.isEmpty())
            {
                String insertBases = breakendInfo.insertedBases;
                if(pairedBreakend.IsRefReversed)
                    insertBases = reverseComplementBases(insertBases);

                BaseSeqViewModel insertSeqViewModel = BaseSeqViewModel.fromStr(insertBases, baseIdx);
                viewRegion = new BaseRegion(baseIdx, baseIdx + breakendInfo.insertedBases.length() - 1);
                BaseRegion assemblyRegion = new BaseRegion(
                        breakendInfo.assemblyRegion.end() + 1, breakendInfo.assemblyRegion.end() + breakendInfo.insertedBases.length());
                baseIdx += breakendInfo.insertedBases.length();
                refViewModelPair.add(new SegmentViewModel(null, null, null, null, viewRegion, null, null, false, false, insertSeqViewModel, true, 0, assemblyRegion));
                boxWidth += breakendInfo.insertedBases.length() + 2 * BOX_PADDING;
            }
        }

        return new PairedRefViewModel(refViewModelPair, baseIdx, boxWidth);
    }

    private static class PairedBreakendInfo
    {
        public final Breakend Breakend;
        public final boolean IsRefReversed;
        public final boolean IsAssemblyReversed;
        public final int AssemblySeqStart;
        public final int AssemblySeqEnd;

        public PairedBreakendInfo(final Breakend breakend)
        {
            Breakend = breakend;
            IsRefReversed = breakend.Orient == REVERSE;
            IsAssemblyReversed = breakend.segments().get(0).Alignment.orientation() == REVERSE;
            AssemblySeqStart = breakend.segments().get(0).Alignment.sequenceStart();
            AssemblySeqEnd = breakend.segments().get(0).Alignment.sequenceEnd();
        }

        @Override
        public boolean equals(final Object o)
        {
            if(!(o instanceof final PairedBreakendInfo that)) return false;
            return AssemblySeqStart == that.AssemblySeqStart && AssemblySeqEnd == that.AssemblySeqEnd;
        }

        @Override
        public int hashCode()
        {
            return AssemblySeqStart + 31 * AssemblySeqEnd;
        }
    }

    private static List<PairedSegmentViewModel> getRefViewModel(final AssemblyConfig config, final AssemblyAlignment assemblyAlignment)
    {
        List<Breakend> breakends = assemblyAlignment.breakends();

        // TODO(mkcmkc): remove
        if(breakends.stream().anyMatch(Breakend::isSingle))
            throw new RuntimeException("Got a single breakend");

        NavigableMap<Integer, Breakend> breakendLookup = Maps.newTreeMap();
        for(Breakend breakend : breakends)
        {
            int id = breakend.id();
            if(breakendLookup.containsKey(id))
            {
                SV_LOGGER.error("Repeated breakend id encountered");
                System.exit(1);
            }

            breakendLookup.put(id, breakend);
        }

        Queue<Breakend> breakedQueue = new PriorityQueue<>(
                Comparator.comparingInt((Breakend x) -> x.segments().get(0).Alignment.sequenceStart())
                        .thenComparingInt((Breakend x) -> x.segments().get(0).Alignment.sequenceEnd()));
        breakedQueue.addAll(breakends);

        List<Pair<PairedBreakendInfo, PairedBreakendInfo>> breakendPairs = Lists.newArrayList();
        Set<Integer> seenIds = Sets.newHashSet();
        while(!breakedQueue.isEmpty())
        {
            Breakend breakend = breakedQueue.poll();
            if(seenIds.contains(breakend.id()))
                continue;

            seenIds.add(breakend.id());
            Breakend mateBreakend = breakend.otherBreakend();
            seenIds.add(mateBreakend.id());

            PairedBreakendInfo pairedBreakendInfo = new PairedBreakendInfo(breakend);
            PairedBreakendInfo matePairedBreakendInfo = new PairedBreakendInfo(mateBreakend);
            breakendPairs.add(Pair.of(pairedBreakendInfo, matePairedBreakendInfo));
        }

        String fullAssemblySeq = assemblyAlignment.fullSequence();
        List<PairedSegmentViewModel> refViewModel = Lists.newArrayList();
        int baseIdx = 0;
        int prevBoxWidth = 0;
        for(Pair<PairedBreakendInfo, PairedBreakendInfo> breakendPair : breakendPairs)
        {
            PairedRefViewModel pairedRefViewModel = getPairedRefViewModel(config, baseIdx, fullAssemblySeq, breakendPair);

            Breakend firstBreakend = breakendPair.getLeft().Breakend;
            int readStartOffset = baseIdx - firstBreakend.segments().get(0).Alignment.sequenceStart();
            refViewModel.add(new PairedSegmentViewModel(pairedRefViewModel.refViewModelPair, readStartOffset, prevBoxWidth));
            prevBoxWidth += pairedRefViewModel.boxWidth;
            baseIdx = pairedRefViewModel.baseIdx;
        }

        return refViewModel;
    }

    private DomContent renderRef()
    {
        int totalBoxWidth = 0;
        for(PairedSegmentViewModel pairedRefEl : mRefViewModel)
        {
            for(SegmentViewModel refEl : pairedRefEl.viewModels)
                totalBoxWidth += refEl.viewRegion.baseLength() + 2 * BOX_PADDING;
        }

        SVGGraphics2D svgCanvas = new SVGGraphics2D(READ_HEIGHT_PX * totalBoxWidth, READ_HEIGHT_PX);
        AffineTransform initTransform = svgCanvas.getTransform();
        double xBoxOffset = 0.0d;
        for(PairedSegmentViewModel pairedRefEl : mRefViewModel)
        {
            for(SegmentViewModel refEl : pairedRefEl.viewModels)
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
        }

        return rawHtml(svgCanvas.getSVGElement());
    }

    private DomContent renderAssembly()
    {
        int totalBoxWidth = 0;
        for(PairedSegmentViewModel pairedRefEl : mRefViewModel)
        {
            for(SegmentViewModel refEl : pairedRefEl.viewModels)
                totalBoxWidth += refEl.viewRegion.baseLength() + 2 * BOX_PADDING;
        }

        SVGGraphics2D svgCanvas = new SVGGraphics2D(READ_HEIGHT_PX * totalBoxWidth, READ_HEIGHT_PX);
        AffineTransform initTransform = svgCanvas.getTransform();
        double xBoxOffset = 0.0d;
        for(PairedSegmentViewModel pairedRefEl : mRefViewModel)
        {
            for(SegmentViewModel refEl : pairedRefEl.viewModels)
            {
                BaseSeqViewModel viewModel = refEl.assemblyViewModel;
                BaseRegion viewRegion = refEl.viewRegion;
                svgCanvas.setTransform(initTransform);
                renderBaseSeq(svgCanvas, new Point2D.Double(xBoxOffset, 0.0d), READ_HEIGHT_PX, viewRegion, viewModel, false, Maps.newHashMap(), refEl.refViewModel());
                xBoxOffset += viewRegion.baseLength() + 2 * BOX_PADDING;
            }
        }

        return rawHtml(svgCanvas.getSVGElement());
    }

    private DomContent renderCoords()
    {
        double scalingFactor = READ_HEIGHT_PX / BASE_BOX_SIZE;
        double charWidth = getStringBounds(COORD_FONT, "9").getWidth();
        double maxStringWidth = 0.0d;
        int totalBoxWidth = 0;
        for(PairedSegmentViewModel pairedRefEl : mRefViewModel)
        {
            for(SegmentViewModel refEl : pairedRefEl.viewModels)
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
        }

        SVGGraphics2D svgCanvas = new SVGGraphics2D(
                scalingFactor * BASE_BOX_SIZE * totalBoxWidth,
                scalingFactor * (maxStringWidth + 2 * charWidth));
        AffineTransform initTransform = svgCanvas.getTransform();
        double xBoxOffset = 0.0d;
        for(PairedSegmentViewModel pairedRefEl : mRefViewModel)
        {
            for(SegmentViewModel refEl : pairedRefEl.viewModels)
            {
                if(refEl.isInsert)
                {
                    xBoxOffset += refEl.viewRegion.baseLength() + 2 * BOX_PADDING;
                    continue;
                }

                BaseRegion refViewRegion = refEl.refViewRegion;
                int centerPosition = refEl.position;
                Point2D.Double canvasSize = new Point2D.Double(
                        scalingFactor * BASE_BOX_SIZE * (refViewRegion.baseLength() + 2 * BOX_PADDING), svgCanvas.getHeight());
                svgCanvas.setTransform(initTransform);
                SvgRender.renderCoords(svgCanvas, new Point2D.Double(xBoxOffset, 0.0d), canvasSize, READ_HEIGHT_PX, refViewRegion, centerPosition, DISPLAY_EVERY_NTH_COORD, maxStringWidth, refEl.isRefReversed);
                xBoxOffset += refViewRegion.baseLength() + 2 * BOX_PADDING;
            }
        }

        return rawHtml(svgCanvas.getSVGElement());
    }

    private DomContent renderChomosomeLabels()
    {
        List<ChrLabel> labels = Lists.newArrayList();
        int xBoxOffset = 0;
        for(PairedSegmentViewModel pairedRefEl : mRefViewModel)
        {
            for(int i = 0; i < pairedRefEl.viewModels.size(); i++)
            {
                SegmentViewModel refEl = pairedRefEl.viewModels.get(i);
                BaseRegion viewRegion = refEl.viewRegion;
                int length = viewRegion.baseLength() + 2 * BOX_PADDING;
                if(refEl.isInsert)
                {
                    labels.add(new ChrLabel(null, 0, new BaseRegion(xBoxOffset, xBoxOffset + length - 1), null, false, null));
                    xBoxOffset += length;
                    continue;
                }

                List<Alignment> alignments = i == 0 ? List.of(RIGHT) : List.of(LEFT);
                labels.add(new ChrLabel(refEl.chromosome, refEl.position, new BaseRegion(xBoxOffset,
                        xBoxOffset + length - 1), refEl.refRegion, refEl.isRefReversed, alignments));
                xBoxOffset += length;
            }
        }

        SVGGraphics2D svgCanvas = SvgRender.renderChrLabels(READ_HEIGHT_PX, labels);
        return rawHtml(svgCanvas.getSVGElement());
    }

    private String getFilename(final String info)
    {
        String sampleId = mConfig.TumorIds.get(0);
        return sampleId + "_" + info.replace(':', '_') + ".html";
    }
}
