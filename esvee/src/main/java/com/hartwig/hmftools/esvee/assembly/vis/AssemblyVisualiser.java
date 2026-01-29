package com.hartwig.hmftools.esvee.assembly.vis;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.codon.Nucleotides.reverseComplementBases;
import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.sv.LineElements.isMobileLineElement;
import static com.hartwig.hmftools.common.sv.SvVcfTags.ALTALN;
import static com.hartwig.hmftools.common.sv.SvVcfTags.ASM_ID;
import static com.hartwig.hmftools.common.sv.SvVcfTags.ASM_INFO;
import static com.hartwig.hmftools.common.sv.SvVcfTags.ASM_LENGTH;
import static com.hartwig.hmftools.common.sv.SvVcfTags.ASM_LINKS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.ASM_SEG_INDEX;
import static com.hartwig.hmftools.common.sv.SvVcfTags.AVG_FRAG_LENGTH;
import static com.hartwig.hmftools.common.sv.SvVcfTags.BE_ASM_ORIENT;
import static com.hartwig.hmftools.common.sv.SvVcfTags.BE_ASM_POS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.BE_ORIENT;
import static com.hartwig.hmftools.common.sv.SvVcfTags.CIPOS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.DISC_FRAGS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.HOMSEQ;
import static com.hartwig.hmftools.common.sv.SvVcfTags.IHOMPOS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.INSALN;
import static com.hartwig.hmftools.common.sv.SvVcfTags.LINE_SITE;
import static com.hartwig.hmftools.common.sv.SvVcfTags.MATE_ID;
import static com.hartwig.hmftools.common.sv.SvVcfTags.MAX_LOCAL_REPEAT;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SEG_ALIGN_LENGTH;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SEG_ID;
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
import static com.hartwig.hmftools.esvee.assembly.alignment.AlternativeAlignment.toVcfTag;
import static com.hartwig.hmftools.esvee.assembly.output.VcfWriter.buildAlleleInfo;
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
import static j2html.TagCreator.table;
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
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Queue;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.vis.BaseSeqViewModel;
import com.hartwig.hmftools.common.vis.CssBuilder;
import com.hartwig.hmftools.common.vis.CssSize;
import com.hartwig.hmftools.common.vis.SvgRender;
import com.hartwig.hmftools.common.vis.SvgRender.ChrLabel_;
import com.hartwig.hmftools.common.vis.SvgUtil;
import com.hartwig.hmftools.common.vis.SvgUtil.Alignment;
import com.hartwig.hmftools.esvee.assembly.AssemblyConfig;
import com.hartwig.hmftools.esvee.assembly.alignment.AlignData;
import com.hartwig.hmftools.esvee.assembly.alignment.AlternativeAlignment;
import com.hartwig.hmftools.esvee.assembly.alignment.AssemblyAlignment;
import com.hartwig.hmftools.esvee.assembly.alignment.Breakend;
import com.hartwig.hmftools.esvee.assembly.alignment.BreakendSegment;
import com.hartwig.hmftools.esvee.assembly.alignment.BreakendSupport;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;

import org.apache.commons.lang3.NotImplementedException;
import org.jetbrains.annotations.Nullable;
import org.jfree.svg.SVGGraphics2D;

import htsjdk.samtools.CigarElement;
import htsjdk.variant.variantcontext.Genotype;
import j2html.tags.DomContent;
import j2html.tags.specialized.TdTag;

public class AssemblyVisualiser
{
    private static final String MISSING_FIELD = "missing";

    public record SegmentViewModel_(String chromosome, @Nullable Integer position, BaseRegion refRegion, BaseRegion viewRegion,
                                   BaseRegion refViewRegion, BaseSeqViewModel refViewModel, boolean isRefReversed,
				   BaseSeqViewModel assemblyViewModel, boolean isInsert, int leftDelLength) {}

    private final AssemblyConfig mConfig;
    private final AssemblyAlignment mAssemblyAlignment;
    private final List<SegmentViewModel_> mRefViewModel_;
    private final List<String> mSampleNames;
    private final Map<String, Integer> mSampleNameIndex;
    private final Set<String> mTumorIds;

    public AssemblyVisualiser(final AssemblyConfig config, final AssemblyAlignment assemblyAlignment)
    {
        mConfig = config;
        mAssemblyAlignment = assemblyAlignment;
        mRefViewModel_ = getRefViewModel(config, assemblyAlignment);
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

        // TODO: remove
//        for(Breakend breakend : mAssemblyAlignment.breakends())
//        {
//            bodyElements.add(renderBreakendInfo(breakend));
//            bodyElements.add(verticalSpacer);
//        }

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

        // TODO: do not need StringJoiner anymore.
        StringJoiner infoJoiner = new StringJoiner(" ");
        infoJoiner.add(format("%s[normal: %s]", mConfig.TumorIds.get(0), mConfig.ReferenceIds.get(0)));

        // TODO: remove
//        infoJoiner.add(mAssemblyAlignment.toString());

        return div(infoJoiner.toString()).withStyle(strongStyle.toString());
    }

    // TODO: remove
    private DomContent renderBreakendInfo(final Breakend breakend)
    {
        LinkedHashMap<String, String> info = new LinkedHashMap<>();

        List<Genotype> genotypes = Lists.newArrayList();
        int totalSplitFrags = 0;
        int totalDiscFrags = 0;
        for(String sampleId : mSampleNames)
        {
            int sampleSupportIndex = mSampleNameIndex.get(sampleId);
            BreakendSupport breakendSupport = breakend.sampleSupport().get(sampleSupportIndex);
            genotypes.add(buildGenotype(sampleId, breakendSupport));

            totalSplitFrags += breakendSupport.SplitFragments;
            totalDiscFrags += breakendSupport.DiscordantFragments;
        }

        info.put("Breakend", breakend.toString());
        info.put("Id", String.valueOf(breakend.id()));
        info.put("Chromosome", breakend.Chromosome);
        info.put("Position", String.valueOf(breakend.Position));
        info.put("Qual", String.valueOf(breakend.calcSvQual()));
        info.put("Alleles", String.valueOf(buildAlleleInfo(mConfig, breakend)));
        info.put("Genotypes", String.valueOf(genotypes));

        if(!breakend.isSingle())
        {
            int otherBreakendId = breakend.otherBreakend().id();
            info.put(MATE_ID, String.valueOf(otherBreakendId));

            // always write lower first
            int lowerBreakendId = otherBreakendId > breakend.id() ? breakend.id() : otherBreakendId;
            int upperBreakendId = lowerBreakendId == breakend.id() ? otherBreakendId : breakend.id();
            String svId = format("%d_%d", lowerBreakendId, upperBreakendId);
            info.put(SV_ID, svId);
        }

        info.put(SV_TYPE, breakend.svType().name());

        if(breakend.Homology.exists())
        {
            info.put(CIPOS, format("%d,%d", breakend.Homology.ExactStart, breakend.Homology.ExactEnd));
            info.put(IHOMPOS, format("%d,%d", breakend.Homology.InexactStart, breakend.Homology.InexactEnd));
            if(!breakend.Homology.Homology.isEmpty())
                info.put(HOMSEQ, breakend.Homology.Homology);
        }
        else
        {
            info.put(CIPOS, "0,0");
            info.put(IHOMPOS, "0,0");
        }

        info.put(SPLIT_FRAGS, String.valueOf(totalSplitFrags));
        info.put(DISC_FRAGS, String.valueOf(totalDiscFrags));
        info.put(TOTAL_FRAGS, String.valueOf(totalSplitFrags + totalDiscFrags));

        // for SGLs with too few fragments with a valid length, set the calculated value to zero
        if(breakend.isSingle() && breakend.validFragmentLengthPercent() < BREAKEND_REQ_VALID_FRAGMENT_LENGTH_PERC)
            info.put(AVG_FRAG_LENGTH, "0");
        else
            info.put(AVG_FRAG_LENGTH, String.valueOf(breakend.averageFragmentLength()));

        List<AlternativeAlignment> altAlignments = breakend.alternativeAlignments();
        if(!altAlignments.isEmpty())
            info.put(INSALN, toVcfTag(altAlignments));

        List<AlternativeAlignment> lowQualAltAlignments = breakend.lowQualAltAlignments();
        if(!lowQualAltAlignments.isEmpty())
            info.put(ALTALN, toVcfTag(lowQualAltAlignments));

        AssemblyAlignment assemblyAlignment = breakend.assembly();

        info.put(ASM_ID, String.valueOf(assemblyAlignment.id()));
        info.put(ASM_INFO, assemblyAlignment.info(10));
        info.put(ASM_LENGTH, String.valueOf(assemblyAlignment.fullSequenceLength()));

        if(assemblyAlignment.assemblies().stream().anyMatch(x -> x.hasLineSequence())
                && isMobileLineElement(breakend.Orient, breakend.InsertedBases))
            info.put(LINE_SITE, "TRUE");
        else
            info.put(LINE_SITE, "FALSE");


        List<BreakendSegment> segments = breakend.segments();

        info.put(ASM_SEG_INDEX, segments.stream().map(x -> String.valueOf(x.Index)).collect(Collectors.joining(VCF_ITEM_DELIM)));
        info.put(BE_ASM_POS, segments.stream().map(x -> String.valueOf(x.SequenceIndex)).collect(Collectors.joining(VCF_ITEM_DELIM)));
        info.put(BE_ORIENT, breakend.Orient.name());
        info.put(BE_ASM_ORIENT, segments.stream().map(x -> x.Alignment.orientation().name()).collect(Collectors.joining(VCF_ITEM_DELIM)));
        info.put(SEG_ID, segments.stream().map(x -> x.uniqueId()).collect(Collectors.joining(VCF_ITEM_DELIM)));

        // NOTE: this is used by Linx to form assembly TIs
        if(!breakend.facingBreakends().isEmpty())
            info.put(ASM_LINKS, breakend.facingBreakends().stream().map(x -> String.valueOf(x.id())).collect(Collectors.joining(VCF_ITEM_DELIM)));

        info.put(SEG_ALIGN_LENGTH, segments.stream().map(x -> String.valueOf(x.Alignment.alignedBases())).collect(Collectors.joining(VCF_ITEM_DELIM)));
        info.put(SEG_MAPQ, String.valueOf(segments.stream().mapToInt(x -> x.Alignment.mapQual()).max().orElse(0)));
        info.put(SEG_SCORE, String.valueOf(segments.stream().mapToInt(x -> x.Alignment.score()).max().orElse(0)));
        info.put(SEG_REPEAT_LENGTH, String.valueOf(segments.stream().mapToInt(x -> x.Alignment.adjustedAlignment()).max().orElse(0)));
        info.put(UNIQUE_FRAG_POSITIONS, String.valueOf(breakend.uniqueFragmentPositionCount()));
        info.put(THREE_PRIME_RANGE, Arrays.stream(breakend.readOrientationRange()).mapToObj(String::valueOf).collect(Collectors.joining(",")));

        if(breakend.maxLocalRepeat() > 0)
            info.put(MAX_LOCAL_REPEAT, String.valueOf(breakend.maxLocalRepeat()));

        CssBuilder readTableStyle = CssBuilder.EMPTY.borderSpacing(CssSize.ZERO);
        List<DomContent> infoRows = Lists.newArrayList();
        for(Map.Entry<String, String> entry : info.entrySet())
            infoRows.add(tr(td(entry.getKey() + ":"), td(entry.getValue())));

        return div(styledTable(infoRows, readTableStyle));
    }

    private static void renderInfoTableSection(final List<DomContent> rows, final String sectionHeader, final LinkedHashMap<String, List<String>> info, final CssBuilder leftHeaderStyle, final CssBuilder firstRowCellStyle, final CssBuilder cellStyle)
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
        CssBuilder firstRowCellStyle = CssBuilder.EMPTY.paddingLeft(CssSize.em(0.5)).paddingRight(CssSize.em(0.5)).merge(firstRowBorderStyle);
        CssBuilder leftHeaderStyle = CssBuilder.EMPTY.backgroundColor(Color.LIGHT_GRAY).writingMode("vertical-rl").transform("rotate(180deg)").fontWeight("bold").textAlign("center").merge(firstRowBorderStyle);

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
                sampleSupportInfo.computeIfAbsent(TOTAL_FRAGS, k -> Lists.newArrayList()).add(String.valueOf(totalDiscordantFrags + totalSplitFrags));
                sampleSupportInfo.computeIfAbsent(DISC_FRAGS, k -> Lists.newArrayList()).add(String.valueOf(totalDiscordantFrags));
                sampleSupportInfo.computeIfAbsent("DP", k -> Lists.newArrayList()).add(String.valueOf(tumorGenotype.getDP()));

                // TODO: other fields
                sampleSupportInfo.computeIfAbsent("REF", k -> Lists.newArrayList()).add("TODO");
                sampleSupportInfo.computeIfAbsent("REFPAIR", k -> Lists.newArrayList()).add("TODO");

                sampleSupportInfo.computeIfAbsent("SB", k -> Lists.newArrayList()).add(format("%.3f", tumorGenotype.getExtendedAttribute("SB")));
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
        CssBuilder firstRowCellStyle = CssBuilder.EMPTY.paddingLeft(CssSize.em(0.5)).paddingRight(CssSize.em(0.5)).merge(firstRowBorderStyle);
        CssBuilder leftHeaderStyle = CssBuilder.EMPTY.backgroundColor(Color.LIGHT_GRAY).writingMode("vertical-rl").transform("rotate(180deg)").fontWeight("bold").textAlign("center").merge(firstRowBorderStyle);

        LinkedHashMap<String, List<String>> idInfo = new LinkedHashMap<>();
        idInfo.put("Id", breakends.stream().map(x -> String.valueOf(x.id())).toList());

        for(Breakend breakend : breakends)
        {
            int otherBreakendId = breakend.otherBreakend().id();
            idInfo.computeIfAbsent(MATE_ID, k -> Lists.newArrayList()).add(String.valueOf(otherBreakendId));

            // always write lower first
            int lowerBreakendId = otherBreakendId > breakend.id() ? breakend.id() : otherBreakendId;
            int upperBreakendId = lowerBreakendId == breakend.id() ? otherBreakendId : breakend.id();
            String svId = format("%d_%d", lowerBreakendId, upperBreakendId);
            idInfo.computeIfAbsent(SV_ID, k -> Lists.newArrayList()).add(svId);
        }

        LinkedHashMap<String, List<String>> segmentPropertiesInfo = new LinkedHashMap<>();
        for(Breakend breakend : breakends)
        {
            List<BreakendSegment> segments = breakend.segments();
            segmentPropertiesInfo.computeIfAbsent(BE_ASM_POS, k -> Lists.newArrayList()).add(
                    segments.stream().map(x -> String.valueOf(x.SequenceIndex)).collect(Collectors.joining(VCF_ITEM_DELIM)));

            segmentPropertiesInfo.computeIfAbsent(SEG_ALIGN_LENGTH, k -> Lists.newArrayList()).add(segments.stream().map(x -> String.valueOf(x.Alignment.alignedBases())).collect(Collectors.joining(VCF_ITEM_DELIM)));
            segmentPropertiesInfo.computeIfAbsent(SEG_MAPQ, k -> Lists.newArrayList()).add(String.valueOf(segments.stream().mapToInt(x -> x.Alignment.mapQual()).max().orElse(0)));
            segmentPropertiesInfo.computeIfAbsent(SEG_SCORE, k -> Lists.newArrayList()).add(String.valueOf(segments.stream().mapToInt(x -> x.Alignment.score()).max().orElse(0)));
            segmentPropertiesInfo.computeIfAbsent(SEG_REPEAT_LENGTH, k -> Lists.newArrayList()).add(String.valueOf(segments.stream().mapToInt(x -> x.Alignment.adjustedAlignment()).max().orElse(0)));
        }

        LinkedHashMap<String, List<String>> siteQualInfo = new LinkedHashMap<>();
        for(Breakend breakend : breakends)
        {
            if(breakend.isSingle() && breakend.validFragmentLengthPercent() < BREAKEND_REQ_VALID_FRAGMENT_LENGTH_PERC)
                siteQualInfo.computeIfAbsent(AVG_FRAG_LENGTH, k -> Lists.newArrayList()).add("0");
            else
                siteQualInfo.computeIfAbsent(AVG_FRAG_LENGTH, k -> Lists.newArrayList()).add(String.valueOf(breakend.averageFragmentLength()));

            siteQualInfo.computeIfAbsent(UNIQUE_FRAG_POSITIONS, k -> Lists.newArrayList()).add(String.valueOf(breakend.uniqueFragmentPositionCount()));
            siteQualInfo.computeIfAbsent(THREE_PRIME_RANGE, k -> Lists.newArrayList()).add(Arrays.stream(breakend.readOrientationRange()).mapToObj(String::valueOf).collect(Collectors.joining(",")));
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
                    ReadViewModel readViewModel = ReadViewModel.create(mRefViewModel_, read, assembly);
                    if(readViewModel == null)
                        continue;

                    readViewModels.add(readViewModel);
                }
                catch(Exception e)
                {
                    e.printStackTrace();
                    throw e;
                }
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

    private record BreakendInfo(String chromosome, int pos, Orientation orient, List<CigarElement> cigarElements, String insertedBases,
                                BaseRegion refRegion, String assemblySeq, String refSeq, boolean isRefReversed, int alignmentSequenceEnd,
                                int leftDelLength) {}

    private static BreakendInfo extractBreakendInfo(
            final RefGenomeInterface refGenome, final Breakend breakend, final String fullAssemblySeq, final boolean isRefReversed,
            @Nullable final Integer lastSeqEnd)
    {
        AlignData alignment = breakend.segments().get(0).Alignment;

        // TODO: remove
        System.out.println("*** alignment range = " + alignment.sequenceStart() + "-" + alignment.sequenceEnd());

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

        if(isRefReversed)
            Collections.reverse(cigarEls);

        String chromosome = alignment.refLocation().chromosome();
        BaseRegion refRegion = alignment.refLocation().baseRegion();
        int position = breakend.Orient == FORWARD ? refRegion.end() : refRegion.start();
        String refSeq = refGenome.getBaseString(chromosome, refRegion.start(), refRegion.end());
        if(isRefReversed)
            refSeq = reverseComplementBases(refSeq);

        return new BreakendInfo(chromosome, position, breakend.Orient, cigarEls, breakend.InsertedBases, refRegion, assemblySeq, refSeq, isRefReversed, alignment.sequenceEnd(), leftDelLength);
    }

    private static List<SegmentViewModel_> getRefViewModel(final AssemblyConfig config, final AssemblyAlignment assemblyAlignment)
    {
        // TODO(mkcmkc): remove
        if(assemblyAlignment.breakends().size() < 2)
        {
            throw new RuntimeException("Less than two breakends.");
            //            SV_LOGGER.error("Expected exactly two breakends.");
            //            System.exit(1);
        }

        // TODO(mkcmkc): remove
//        if(assemblyAlignment.breakends().size() >= 4)
//        {
//            throw new RuntimeException("More than 3 breakends.");
//            //            SV_LOGGER.error("Expected exactly two breakends.");
//            //            System.exit(1);
//        }

        // TODO(mkcmkc): remove
        if(assemblyAlignment.breakends().size() == 2)
        {
            throw new RuntimeException("Exactly 2 breakends.");
            //            SV_LOGGER.error("Expected exactly two breakends.");
            //            System.exit(1);
        }


        List<SegmentViewModel_> refViewModel_ = Lists.newArrayList();
        int breakendCount = assemblyAlignment.breakends().size();
        String fullAssemblySeq = assemblyAlignment.fullSequence();
        Integer lastSequenceEnd = null;
        int baseIdx = 0;
        for(int i_ = 0; i_ < breakendCount; i_++)
        {
            Breakend breakend = assemblyAlignment.breakends().get(i_);
            boolean isRefReversed_ = (i_ < breakendCount - 1 && breakend.Orient == REVERSE) || (i_ == breakendCount - 1 && breakend.Orient == FORWARD);
            BreakendInfo breakendInfo = extractBreakendInfo(config.RefGenome, breakend, fullAssemblySeq, isRefReversed_, lastSequenceEnd);
            lastSequenceEnd = breakendInfo.alignmentSequenceEnd;

            // TODO: remove
            System.out.println("*** " + breakend);
            System.out.println("*** " + breakendInfo);

            BaseSeqViewModel refSeqViewModel = BaseSeqViewModel.fromStr(breakendInfo.refSeq, baseIdx);
            BaseSeqViewModel assemblySeqViewModel = BaseSeqViewModel.fromStringWithCigar(
                    breakendInfo.assemblySeq, breakendInfo.cigarElements, baseIdx);

            BaseRegion viewRegion;
            if(i_ == 0)
            {
                int viewRegionEnd = baseIdx + breakendInfo.refRegion.baseLength() - 1;
                int viewRegionStart = max(baseIdx, viewRegionEnd - VIEW_REGION_SIZE + 1);
                viewRegion = new BaseRegion(viewRegionStart, viewRegionEnd);
            }
            else if(i_ == assemblyAlignment.breakends().size() - 1)
            {
                int viewRegionStart = baseIdx;
                int viewRegionEnd = min(baseIdx + breakendInfo.refRegion.baseLength() - 1, viewRegionStart + VIEW_REGION_SIZE - 1);
                viewRegion = new BaseRegion(viewRegionStart, viewRegionEnd);
            }
            else
            {
                viewRegion = new BaseRegion(baseIdx, baseIdx + breakendInfo.refRegion.baseLength() - 1);
            }

            baseIdx += breakendInfo.refRegion.baseLength();

            BaseRegion refViewRegion;
            if(i_ == 0 || i_ == assemblyAlignment.breakends().size() - 1)
            {
                if(breakendInfo.orient == FORWARD)
                {
                    int viewRegionEnd = breakendInfo.refRegion.end();
                    int viewRegionStart = max(breakendInfo.refRegion.start(), viewRegionEnd - VIEW_REGION_SIZE + 1);
                    refViewRegion = new BaseRegion(viewRegionStart, viewRegionEnd);
                }
                else
                {
                    int viewRegionStart = breakendInfo.refRegion.start();
                    int viewRegionEnd = min(breakendInfo.refRegion.end(), viewRegionStart + VIEW_REGION_SIZE - 1);
                    refViewRegion = new BaseRegion(viewRegionStart, viewRegionEnd);
                }
            }
            else
            {
                refViewRegion = breakendInfo.refRegion;
            }

            refViewModel_.add(new SegmentViewModel_(breakendInfo.chromosome, breakendInfo.pos, breakendInfo.refRegion, viewRegion, refViewRegion, refSeqViewModel, isRefReversed_, assemblySeqViewModel, false, breakendInfo.leftDelLength));
            if(i_ < assemblyAlignment.breakends().size() - 1 && breakendInfo.insertedBases != null && !breakendInfo.insertedBases.isEmpty())
            {
                String insertBases = breakendInfo.insertedBases;
                if(isRefReversed_)
                    insertBases = reverseComplementBases(insertBases);

                BaseSeqViewModel insertSeqViewModel = BaseSeqViewModel.fromStr(insertBases, baseIdx);
                viewRegion = new BaseRegion(baseIdx, baseIdx + breakendInfo.insertedBases.length() - 1);
                baseIdx += breakendInfo.insertedBases.length();
                refViewModel_.add(new SegmentViewModel_(null, null, null, viewRegion, null, null, false, insertSeqViewModel, true, 0));
            }
        }

        return refViewModel_;
    }

    private DomContent renderRef()
    {
        int totalBoxWidth = 0;
        for(SegmentViewModel_ refEl : mRefViewModel_)
            totalBoxWidth += refEl.viewRegion.baseLength() + 2 * BOX_PADDING;

        SVGGraphics2D svgCanvas = new SVGGraphics2D(READ_HEIGHT_PX * totalBoxWidth, READ_HEIGHT_PX);
        AffineTransform initTransform = svgCanvas.getTransform();
        double xBoxOffset = 0.0d;
        for(SegmentViewModel_ refEl : mRefViewModel_)
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
        for(SegmentViewModel_ refEl : mRefViewModel_)
            totalBoxWidth += refEl.viewRegion.baseLength() + 2 * BOX_PADDING;

        SVGGraphics2D svgCanvas = new SVGGraphics2D(READ_HEIGHT_PX * totalBoxWidth, READ_HEIGHT_PX);
        AffineTransform initTransform = svgCanvas.getTransform();
        double xBoxOffset = 0.0d;
        for(SegmentViewModel_ refEl : mRefViewModel_)
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
        for(SegmentViewModel_ refEl : mRefViewModel_)
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
        for(SegmentViewModel_ refEl : mRefViewModel_)
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

        return rawHtml(svgCanvas.getSVGElement());
    }

    private DomContent renderChomosomeLabels()
    {
        List<ChrLabel_> labels_ = Lists.newArrayList();

        int xBoxOffset = 0;
        for(int i = 0; i < mRefViewModel_.size(); i++)
        {
            SegmentViewModel_ refEl = mRefViewModel_.get(i);
            BaseRegion viewRegion = refEl.viewRegion;
            int length = viewRegion.baseLength() + 2 * BOX_PADDING;
            if(refEl.isInsert)
            {
                labels_.add(new ChrLabel_(null, 0, new BaseRegion(xBoxOffset, xBoxOffset + length - 1), null, false, null));
                xBoxOffset += length;
                continue;
            }

            List<Alignment> alignments;
            if(i == 0)
                alignments = List.of(RIGHT);
            else if(i == mRefViewModel_.size() - 1)
                alignments = List.of(LEFT);
            else
                alignments = List.of(LEFT, RIGHT);

            labels_.add(new ChrLabel_(refEl.chromosome, refEl.position, new BaseRegion(xBoxOffset, xBoxOffset + length - 1), refEl.refRegion, refEl.isRefReversed, alignments));
            xBoxOffset += length;
        }

        // TODO:
        System.out.println("*** labels = " + labels_);

        SVGGraphics2D svgCanvas = SvgRender.renderChrLabels(READ_HEIGHT_PX, labels_);
        return rawHtml(svgCanvas.getSVGElement());
    }

    private String getFilename(final String info)
    {
        String sampleId = mConfig.TumorIds.get(0);
        return sampleId + "_" + info.replace(':', '_') + ".html";
    }
}
