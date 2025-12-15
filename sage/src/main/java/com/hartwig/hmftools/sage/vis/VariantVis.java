package com.hartwig.hmftools.sage.vis;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;
import static java.lang.String.format;
import static java.util.Map.entry;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome.CHR_PREFIX;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.variant.SageVcfTags.AVG_RECALIBRATED_BASE_QUAL;
import static com.hartwig.hmftools.common.variant.SageVcfTags.MIN_COORDS_COUNT;
import static com.hartwig.hmftools.common.variant.SageVcfTags.UMI_TYPE_COUNTS;
import static com.hartwig.hmftools.common.vis.BaseSeqViewModel.fromStr;
import static com.hartwig.hmftools.common.vis.ColorUtil.DARK_BLUE;
import static com.hartwig.hmftools.common.vis.HtmlUtil.BASE_FONT_STYLE;
import static com.hartwig.hmftools.common.vis.HtmlUtil.JQUERY_SCRIPT;
import static com.hartwig.hmftools.common.vis.HtmlUtil.getJavascript;
import static com.hartwig.hmftools.common.vis.HtmlUtil.renderReadInfoTable;
import static com.hartwig.hmftools.common.vis.HtmlUtil.styledTable;
import static com.hartwig.hmftools.common.vis.SvgRender.renderBaseSeq;
import static com.hartwig.hmftools.common.vis.SvgRender.renderCoords;
import static com.hartwig.hmftools.common.vis.SvgRender.renderGeneData;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.vcf.VcfTags.AVG_READ_MAP_QUALITY;
import static com.hartwig.hmftools.sage.vcf.VcfTags.AVG_SEQ_TECH_BASE_QUAL;
import static com.hartwig.hmftools.sage.vcf.VcfTags.FRAG_STRAND_BIAS;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_STRAND_BIAS;
import static com.hartwig.hmftools.sage.vis.AminoAcidUtil.getAltGeneRegionViewModels;
import static com.hartwig.hmftools.sage.vis.AminoAcidUtil.getGeneRegionLabel;
import static com.hartwig.hmftools.sage.vis.AminoAcidUtil.getGeneRegions;
import static com.hartwig.hmftools.sage.vis.AminoAcidUtil.getRefGeneRegionViewModels;
import static com.hartwig.hmftools.sage.vis.ReadTableColumn.FINAL_BASE_QUAL_COL;
import static com.hartwig.hmftools.sage.vis.ReadTableColumn.FINAL_MAP_QUAL_COL;
import static com.hartwig.hmftools.sage.vis.ReadTableColumn.FINAL_QUAL_COL;
import static com.hartwig.hmftools.sage.vis.ReadTableColumn.MAP_QUAL_COL;
import static com.hartwig.hmftools.sage.vis.ReadTableColumn.MATE_TYPE_COL;
import static com.hartwig.hmftools.sage.vis.ReadTableColumn.ORIENTATION_COL;
import static com.hartwig.hmftools.sage.vis.ReadTableColumn.SEQ_TECH_BASE_QUAL_COL;
import static com.hartwig.hmftools.sage.vis.SageVisConstants.AA_VARIANT_TYPE_IDX;
import static com.hartwig.hmftools.sage.vis.SageVisConstants.DISPLAY_EVERY_NTH_COORD;
import static com.hartwig.hmftools.sage.vis.SageVisConstants.GENE_NAME_IDX;
import static com.hartwig.hmftools.sage.vis.SageVisConstants.HGVS_INDEX;
import static com.hartwig.hmftools.sage.vis.SageVisConstants.IMPACT_KEY;
import static com.hartwig.hmftools.sage.vis.SageVisConstants.MAX_READ_UPPER_LIMIT;
import static com.hartwig.hmftools.sage.vis.SageVisConstants.READ_HEIGHT_PX;
import static com.hartwig.hmftools.sage.vis.SageVisConstants.TRANSCRIPT_NAME_IDX;

import static j2html.TagCreator.body;
import static j2html.TagCreator.div;
import static j2html.TagCreator.header;
import static j2html.TagCreator.html;
import static j2html.TagCreator.rawHtml;
import static j2html.TagCreator.span;
import static j2html.TagCreator.td;
import static j2html.TagCreator.tr;

import java.awt.Color;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.EnumMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.SortedSet;
import java.util.StringJoiner;
import java.util.concurrent.ThreadLocalRandom;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.variant.SimpleVariant;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.common.vis.BaseSeqViewModel;
import com.hartwig.hmftools.common.vis.BaseViewModel;
import com.hartwig.hmftools.common.vis.CssBuilder;
import com.hartwig.hmftools.common.vis.CssSize;
import com.hartwig.hmftools.common.vis.GeneRegionViewModel;
import com.hartwig.hmftools.common.vis.SvgRender;
import com.hartwig.hmftools.sage.ReferenceData;
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.sage.common.ReadContextMatch;
import com.hartwig.hmftools.sage.common.SageVariant;
import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;
import com.hartwig.hmftools.sage.quality.QualityScores;
import com.hartwig.hmftools.sage.sync.FragmentData;
import com.hartwig.hmftools.sage.vis.ReadTableColumn.ContentAndStyle;

import org.jetbrains.annotations.Nullable;
import org.jfree.svg.SVGGraphics2D;

import htsjdk.samtools.SAMRecord;
import htsjdk.variant.variantcontext.VariantContext;
import j2html.tags.DomContent;
import j2html.tags.specialized.TdTag;

public class VariantVis
{
    private static final List<ReadContextMatch> SORTED_MATCH_TYPES = Arrays.stream(ReadContextMatch.values())
            .sorted(Comparator.comparingInt(VariantVis::visSortKey))
            .toList();

    private static int visSortKey(final ReadContextMatch type)
    {
        return switch(type)
        {
            case FULL -> 0;
            case PARTIAL_CORE -> 1;
            case REALIGNED -> 2;
            case CORE -> 3;
            case REF -> 3;
            default -> 5;
        };
    }

    private static final Map<String, SortedSet<String>> VARIANT_INDEXED_BASES_MAP = Maps.newConcurrentMap();

    private final VisConfig mConfig;

    private final String mSample;
    private final SimpleVariant mVariant;
    private final VariantTier mVariantTier;
    private final EnumMap<ReadContextMatch, List<ReadEvidenceRecord>> mReadEvidenceRecordsByType;
    private final VariantReadContext mReadContext;
    private final BaseRegion mViewRegion;
    private final Map<Integer, List<SvgRender.BoxBorder>> mContextBorders;
    private final BaseSeqViewModel mRefViewModel;
    private final BaseSeqViewModel mContextViewModel;
    private final EnumMap<ReadContextMatch, Integer> mReadCountByType;
    private final String mVariantKey;
    private final String mIndexedBasesKey;
    private final RefGenomeSource mRefGenome;
    private final RefGenomeCoordinates mRefGenCoords;

    private int mReadCount;

    public VariantVis(
            final SageConfig config, final String sample, final SimpleVariant variant, final VariantReadContext readContext,
            final VariantTier variantTier)
    {
        mConfig = config.Visualiser;
        mSample = sample;
        mVariant = variant;
        mVariantTier = variantTier;
        mReadEvidenceRecordsByType = Maps.newEnumMap(ReadContextMatch.class);
        mReadContext = readContext;
        mVariantKey = mVariant.chromosome() + "_" + mVariant.position() + "_" + mVariant.ref() + "_" + mVariant.alt();

        int coreStart = mReadContext.CorePositionStart;
        int coreEnd = mReadContext.CorePositionEnd;
        int flankStart = coreStart - mReadContext.leftFlankLength();
        int flankEnd = coreEnd + mReadContext.rightFlankLength();
        mViewRegion = new BaseRegion(coreStart - SageVisConstants.READ_EXTEND_LENGTH, coreEnd + SageVisConstants.READ_EXTEND_LENGTH);

        mContextBorders = Map.ofEntries(
                entry(mVariant.Position,
                        Lists.newArrayList(
                                new SvgRender.BoxBorder(SvgRender.BorderLocation.LEFT, Color.BLACK),
                                new SvgRender.BoxBorder(SvgRender.BorderLocation.RIGHT, Color.BLACK))),
                entry(coreStart,
                        Lists.newArrayList(
                                new SvgRender.BoxBorder(SvgRender.BorderLocation.LEFT, Color.BLUE))),
                entry(coreEnd,
                        Lists.newArrayList(
                                new SvgRender.BoxBorder(SvgRender.BorderLocation.RIGHT, Color.BLUE))),
                entry(flankStart,
                        Lists.newArrayList(
                                new SvgRender.BoxBorder(SvgRender.BorderLocation.LEFT, DARK_BLUE))),
                entry(flankEnd,
                        Lists.newArrayList(
                                new SvgRender.BoxBorder(SvgRender.BorderLocation.RIGHT, DARK_BLUE)))
        );

        mRefGenCoords = config.RefGenVersion == V37 ? RefGenomeCoordinates.COORDS_37 : RefGenomeCoordinates.COORDS_38;
        int chromosomeLength = mRefGenCoords.length(mVariant.chromosome());
        int refPosStart = max(1, mViewRegion.start());
        int refPosEnd = min(chromosomeLength, mViewRegion.end());
        mRefGenome = loadRefGenome(config.RefGenomeFile);
        String refBases = mRefGenome.getBaseString(mVariant.chromosome(), refPosStart, refPosEnd);
        mRefViewModel = fromStr(refBases, refPosStart);

        mContextViewModel = contextViewModelFromVariant(mReadContext, mVariant.ref(), mVariant.alt());

        StringJoiner indexedBasesKeyBuilder = new StringJoiner("_");
        indexedBasesKeyBuilder.add(String.valueOf(mReadContext.VarIndex));
        indexedBasesKeyBuilder.add(String.valueOf(mReadContext.CoreIndexStart));
        indexedBasesKeyBuilder.add(String.valueOf(mReadContext.CoreIndexEnd));
        indexedBasesKeyBuilder.add(String.valueOf(mReadContext.leftFlankLength()));
        indexedBasesKeyBuilder.add(new String(mReadContext.ReadBases));
        mIndexedBasesKey = indexedBasesKeyBuilder.toString();

        mReadCountByType = Maps.newEnumMap(ReadContextMatch.class);
        mReadCount = 0;
        VARIANT_INDEXED_BASES_MAP.computeIfAbsent(mVariantKey, k -> Sets.newTreeSet()).add(mIndexedBasesKey);
    }

    private static BaseSeqViewModel contextViewModelFromVariant(final VariantReadContext readContext, final String ref, final String alt)
    {
        String rawBases = readContext.readBases();
        int posStart = readContext.variant().Position - readContext.VarIndex;
        if(ref.length() == alt.length())
        {
            return fromStr(rawBases, posStart);
        }

        // del
        if(ref.length() > alt.length())
        {
            int delLen = ref.length() - alt.length();

            // alt is single char
            List<BaseViewModel> bases = Lists.newArrayList();
            for(int i = 0; i <= readContext.VarIndex; ++i)
            {
                bases.add(new BaseViewModel(rawBases.charAt(i)));
            }

            for(int i = 0; i < delLen; ++i)
            {
                bases.add(BaseViewModel.createDelBase());
            }

            for(int i = readContext.VarIndex + 1; i < readContext.totalLength(); ++i)
            {
                bases.add(new BaseViewModel(rawBases.charAt(i)));
            }

            return new BaseSeqViewModel(bases, posStart, null, null);
        }

        // ins
        int insLen = alt.length() - ref.length();

        // ref is single char
        List<BaseViewModel> bases = Lists.newArrayList();
        for(int i = 0; i <= readContext.VarIndex; ++i)
        {
            bases.add(new BaseViewModel(rawBases.charAt(i)));
        }

        bases.get(bases.size() - 1).incRightInsertCount(insLen);

        for(int i = readContext.VarIndex + insLen + 1; i < readContext.totalLength(); ++i)
        {
            bases.add(new BaseViewModel(rawBases.charAt(i)));
        }

        return new BaseSeqViewModel(bases, posStart, null, null);
    }

    public static void writeToHtmlFile(final SageVariant sageVariant, final List<String> tumorIds, final List<String> referenceIds,
            final VisConfig config, @Nullable final ReferenceData refData)
    {
        if(config.PassOnly && !sageVariant.isPassing())
            return;

        List<ReadContextCounter> tumorReadCounters = sageVariant.tumorReadCounters();
        List<ReadContextCounter> refReadCounters = sageVariant.referenceReadCounters();

        // can be null if doesn't meet the configured criteria
        List<VariantVis> tumorVis = tumorReadCounters.stream()
                .map(ReadContextCounter::variantVis).filter(Objects::nonNull).toList();

        List<VariantVis> refVis = refReadCounters.stream()
                .map(ReadContextCounter::variantVis).filter(Objects::nonNull).toList();

        if(tumorVis.isEmpty() && refVis.isEmpty())
            return;

        ReadContextCounter firstCounter = !tumorReadCounters.isEmpty() ? tumorReadCounters.get(0) : refReadCounters.get(0);
        VariantVis firstVis = !tumorVis.isEmpty() ? tumorVis.get(0) : refVis.get(0);

        AminoAcidElements aaElements = getAminoAcidsElements(config, refData, firstVis.mViewRegion, sageVariant, firstVis.mRefGenome);

        String aaVariantType = null;
        String geneName = null;
        if(aaElements != null && aaElements.variant != null)
        {
            VariantContext variantContext = aaElements.variant;
            List<String> impact = (List<String>) variantContext.getAttribute(IMPACT_KEY);
            geneName = impact.get(GENE_NAME_IDX);
            aaVariantType = impact.get(AA_VARIANT_TYPE_IDX);
        }

        String sampleId = tumorIds.isEmpty() ? referenceIds.get(0) : tumorIds.get(0);
        String filename = firstVis.getFilename(sampleId, geneName, aaVariantType);

        int tumorTotalReadCount = tumorVis.stream().mapToInt(x -> x.mReadCount).sum();
        int refTotalReadCount = refVis.stream().mapToInt(x -> x.mReadCount).sum();
        int totalReadCount = tumorTotalReadCount + refTotalReadCount;
        if(totalReadCount == 0)
        {
            SG_LOGGER.info("not writing variant vis file {}, because there are no associated reads", filename);
            return;
        }

        tumorVis.forEach(VariantVis::downsampleReadEvidenceRecords);
        refVis.forEach(VariantVis::downsampleReadEvidenceRecords);

        Stream<DomContent> tumorReadTableRows = tumorVis.stream().map(x -> x.renderReads(true, aaElements)).flatMap(Collection::stream);
        Stream<DomContent> refReadTableRows = refVis.stream().map(x -> x.renderReads(false, aaElements)).flatMap(Collection::stream);
        List<DomContent> readTableRows = Stream.concat(tumorReadTableRows, refReadTableRows).collect(Collectors.toList());

        CssBuilder readTableStyle = CssBuilder.EMPTY.borderSpacing(CssSize.ZERO);
        CssBuilder verticalSpacerStyle = CssBuilder.EMPTY.height(CssSize.em(1));

        DomContent readTable = div(styledTable(readTableRows, readTableStyle));
        DomContent verticalSpacer = div().withStyle(verticalSpacerStyle.toString());
        DomContent variantInfo = firstVis.renderVariantInfo(tumorIds, referenceIds);
        DomContent sampleInfo = styledTable(List.of(tr(
                td(renderSampleSummaryAndCountsTable(tumorReadCounters, refReadCounters))
                        .withStyle(CssBuilder.EMPTY.verticalAlign("top")
                                .toString()),
                td(renderSampleQualAndSiteInfo(tumorReadCounters, refReadCounters, firstCounter))
                        .withStyle(CssBuilder.EMPTY.verticalAlign("top")
                                .toString()),
                td(firstVis.renderVariantInfoTable(
                        (int) (-10 * firstCounter.logTqp()),
                        round(round(firstCounter.mapQualFactor() * 10.0d) / 10.0d),
                        sageVariant.nearIndel(),
                        sageVariant.filtersStringSet(),
                        aaElements)).withStyle(CssBuilder.EMPTY.verticalAlign("top").toString()))), CssBuilder.EMPTY);

        String htmlStr = html(
                header(JQUERY_SCRIPT),
                body(
                        variantInfo,
                        verticalSpacer,
                        sampleInfo,
                        readTable,
                        getJavascript()).withStyle(BASE_FONT_STYLE.toString())).render();

        String filePath = (new File(firstVis.mConfig.OutputDir, filename)).toString();

        SG_LOGGER.debug("writing variant vis file: {}", filePath);

        try
        {
            BufferedWriter outputWriter = createBufferedWriter(filePath, false);
            outputWriter.write(htmlStr);
            outputWriter.newLine();
            closeBufferedWriter(outputWriter);
        }
        catch(IOException e)
        {
            SG_LOGGER.error("failed to write output file({}): {}", filePath, e.toString());
            System.exit(1);
        }
    }

    private static List<TdTag> getSampleInfoRowElements(final String key, @Nullable final String tumor, @Nullable final String ref)
    {
        List<TdTag> elements = Lists.newArrayList();
        elements.add(td(key));
        if(tumor != null)
            elements.add(td(tumor));

        if(ref != null)
            elements.add(td(ref));

        return elements;
    }

    private static DomContent renderSampleSummaryAndCountsTable(
            final List<ReadContextCounter> tumorReadCounters, final List<ReadContextCounter> refReadCounters)
    {
        CssBuilder borderStyle = CssBuilder.EMPTY.border(CssSize.px(1.0), "solid", Color.BLACK).borderCollapse("collapse");
        CssBuilder tableStyle = borderStyle.marginRight(CssSize.px(10.0));
        CssBuilder headerStyle = CssBuilder.EMPTY.fontWeight("bold").textAlign("center").backgroundColor(Color.LIGHT_GRAY);
        CssBuilder cellStyle = CssBuilder.EMPTY.paddingLeft(CssSize.em(0.5)).paddingRight(CssSize.em(0.5));

        ReadContextCounter tumorCounter = null;
        int colCount = 1;
        for(ReadContextCounter counter : tumorReadCounters)
        {
            if(counter.variantVis() == null)
                continue;

            tumorCounter = counter;
            colCount++;
            break;
        }

        ReadContextCounter refCounter = null;
        for(ReadContextCounter counter : refReadCounters)
        {
            if(counter.variantVis() == null)
                continue;

            refCounter = counter;
            colCount++;
            break;
        }

        List<List<TdTag>> contentRowsElems = Lists.newArrayList();
        String tumorValue = tumorCounter == null ? null : String.valueOf(tumorCounter.altSupport());
        String refValue = refCounter == null ? null : String.valueOf(refCounter.altSupport());
        List<TdTag> elems = getSampleInfoRowElements("AD", tumorValue, refValue);
        contentRowsElems.add(elems);

        tumorValue = tumorCounter == null ? null : String.valueOf(tumorCounter.depth());
        refValue = refCounter == null ? null : String.valueOf(refCounter.depth());
        elems = getSampleInfoRowElements("DP", tumorValue, refValue);
        contentRowsElems.add(elems);

        tumorValue = tumorCounter == null ? null : format("%.3f", tumorCounter.vaf());
        refValue = refCounter == null ? null : format("%.3f", refCounter.vaf());
        elems = getSampleInfoRowElements("AF", tumorValue, refValue);
        contentRowsElems.add(elems);

        VariantVis tumorVariantVis = tumorCounter == null ? null : tumorCounter.variantVis();
        VariantVis refVariantVis = refCounter == null ? null : refCounter.variantVis();
        for(ReadContextMatch matchType : SORTED_MATCH_TYPES)
        {
            String key = matchType.name();
            tumorValue = tumorVariantVis == null ? null : String.valueOf(tumorVariantVis.mReadCountByType.getOrDefault(matchType, 0));
            refValue = refVariantVis == null ? null : String.valueOf(refVariantVis.mReadCountByType.getOrDefault(matchType, 0));
            elems = getSampleInfoRowElements(key, tumorValue, refValue);
            contentRowsElems.add(elems);
        }

        List<DomContent> rows = Lists.newArrayList();

        tumorValue = tumorCounter == null ? null : "Tumor";
        refValue = refCounter == null ? null : "Normal";
        elems = getSampleInfoRowElements("Info", tumorValue, refValue).stream()
                .map(x -> x.withStyle(cellStyle.merge(borderStyle).toString()))
                .toList();

        rows.add(tr().with(elems).withStyle(headerStyle.toString()));
        for(int i = 0; i < contentRowsElems.size(); i++)
        {
            // spacer row
            if(i == 3)
                rows.add(tr(td("Support By Type").withStyle(cellStyle.noBorder().toString())
                        .attr("colspan", colCount)).withStyle(headerStyle.toString()));

            List<TdTag> contentRowElems = contentRowsElems.get(i);
            contentRowElems.set(0, contentRowElems.get(0).withStyle(cellStyle.merge(borderStyle).toString()));
            for(int j = 1; j < contentRowElems.size(); j++)
                contentRowElems.set(j, contentRowElems.get(j)
                        .withStyle(CssBuilder.EMPTY.textAlign("right").merge(cellStyle.merge(borderStyle)).toString()));

            rows.add(tr().with(contentRowElems));
        }

        DomContent table = styledTable(rows, tableStyle);
        return div(table);
    }

    private static DomContent renderSampleQualAndSiteInfo(final List<ReadContextCounter> tumorReadCounters,
            final List<ReadContextCounter> refReadCounters, final ReadContextCounter firstCounter)
    {
        CssBuilder borderStyle = CssBuilder.EMPTY.border(CssSize.px(1.0), "solid", Color.BLACK).borderCollapse("collapse");
        CssBuilder tableStyle = borderStyle.marginRight(CssSize.px(10.0));
        CssBuilder headerStyle = CssBuilder.EMPTY.fontWeight("bold").textAlign("center").backgroundColor(Color.LIGHT_GRAY);
        CssBuilder cellStyle = CssBuilder.EMPTY.paddingLeft(CssSize.em(0.5)).paddingRight(CssSize.em(0.5)).merge(borderStyle);

        ReadContextCounter tumorCounter = null;
        for(ReadContextCounter counter : tumorReadCounters)
        {
            if(counter.variantVis() == null)
                continue;

            tumorCounter = counter;
            break;
        }

        ReadContextCounter refCounter = null;
        for(ReadContextCounter counter : refReadCounters)
        {
            if(counter.variantVis() == null)
                continue;

            refCounter = counter;
            break;
        }

        double nonAvgEdgeDist = firstCounter.readEdgeDistance().avgDistanceFromEdge();
        double altAvgEdgeDist = firstCounter.readEdgeDistance().avgAltDistanceFromEdge();

        Map<String, List<String>> values = Maps.newLinkedHashMap();

        String tumorValue = tumorCounter == null ? null : String.valueOf(tumorCounter.tumorQuality());
        String refValue = refCounter == null ? null : String.valueOf(refCounter.tumorQuality());
        values.put("RAW_QUAL", List.of(tumorValue, refValue));

        tumorValue = tumorCounter == null ? null : String.valueOf((int) tumorCounter.averageAltRecalibratedBaseQuality());
        refValue = refCounter == null ? null : String.valueOf((int) refCounter.averageAltRecalibratedBaseQuality());
        values.put(AVG_RECALIBRATED_BASE_QUAL, List.of(tumorValue, refValue));

        tumorValue = tumorCounter == null ? null : String.valueOf((int) tumorCounter.averageAltSeqTechBaseQuality());
        refValue = refCounter == null ? null : String.valueOf((int) refCounter.averageAltSeqTechBaseQuality());
        values.put(AVG_SEQ_TECH_BASE_QUAL, List.of(tumorValue, refValue));

        tumorValue = tumorCounter == null
                ? null
                : String.valueOf(tumorCounter.altSupport() > 0
                        ? (int) round(tumorCounter.altMapQualityTotal() / (double) tumorCounter.altSupport())
                        : 0);
        refValue = refCounter == null
                ? null
                : String.valueOf(
                        refCounter.altSupport() > 0 ? (int) round(refCounter.altMapQualityTotal() / (double) refCounter.altSupport()) : 0);
        values.put(AVG_READ_MAP_QUALITY, List.of(tumorValue, refValue));

        tumorValue = tumorCounter == null ? null : format("%.2f", tumorCounter.fragmentStrandBiasAlt().bias());
        refValue = refCounter == null ? null : format("%.2f", refCounter.fragmentStrandBiasAlt().bias());
        values.put(FRAG_STRAND_BIAS, List.of(tumorValue, refValue));

        tumorValue = tumorCounter == null ? null : format("%.2f", tumorCounter.readStrandBiasAlt().bias());
        refValue = refCounter == null ? null : format("%.2f", refCounter.readStrandBiasAlt().bias());
        values.put(READ_STRAND_BIAS, List.of(tumorValue, refValue));

        tumorValue = tumorCounter == null ? null : format("%d-%d", tumorCounter.jitter().shortened(), tumorCounter.jitter().lengthened());
        refValue = refCounter == null ? null : format("%d-%d", refCounter.jitter().shortened(), refCounter.jitter().lengthened());
        values.put("JIT", List.of(tumorValue, refValue));

        tumorValue = tumorCounter == null ? null : String.valueOf(tumorCounter.fragmentCoords().minCount());
        refValue = refCounter == null ? null : String.valueOf(refCounter.fragmentCoords().minCount());
        values.put(MIN_COORDS_COUNT, List.of(tumorValue, refValue));

        tumorValue = tumorCounter == null ? null : Arrays.toString(tumorCounter.consensusTypeCounts()).replace("[", "").replace("]", "");
        refValue = refCounter == null ? null : Arrays.toString(refCounter.consensusTypeCounts()).replace("[", "").replace("]", "");
        values.put(UMI_TYPE_COUNTS, List.of(tumorValue, refValue));

        tumorValue = tumorCounter == null ? null : format("%.2f", nonAvgEdgeDist);
        refValue = refCounter == null ? null : format("%.2f", altAvgEdgeDist);
        values.put("AED", List.of(tumorValue, refValue));

        List<List<TdTag>> contentRowsElems = Lists.newArrayList();
        for(Map.Entry<String, List<String>> entry : values.entrySet())
            contentRowsElems.add(getSampleInfoRowElements(entry.getKey(), entry.getValue().get(0), entry.getValue().get(1)));

        List<DomContent> rows = Lists.newArrayList();
        tumorValue = tumorCounter == null ? null : "Tumor";
        refValue = refCounter == null ? null : "Normal";
        List<TdTag> elems = getSampleInfoRowElements("Info", tumorValue, refValue).stream()
                .map(x -> x.withStyle(cellStyle.toString()))
                .toList();

        rows.add(tr().with(elems).withStyle(headerStyle.toString()));
        for(List<TdTag> contentRowElems : contentRowsElems)
        {
            contentRowElems.set(0, contentRowElems.get(0).withStyle(cellStyle.toString()));
            for(int j = 1; j < contentRowElems.size(); j++)
                contentRowElems.set(j, contentRowElems.get(j).withStyle(cellStyle.textAlign("right").toString()));

            rows.add(tr().with(contentRowElems));
        }

        DomContent table = styledTable(rows, tableStyle);
        return div(table);
    }

    private DomContent renderVariantInfoTable(int totalTumorQuality, double mapQualFactor, boolean nearbyIndel, final Set<String> filters,
            @Nullable final AminoAcidElements aaElements)
    {
        CssBuilder borderStyle = CssBuilder.EMPTY.border(CssSize.px(1.0), "solid", Color.BLACK).borderCollapse("collapse");
        CssBuilder tableStyle = borderStyle;
        CssBuilder headerStyle = CssBuilder.EMPTY.fontWeight("bold").textAlign("center").backgroundColor(Color.LIGHT_GRAY);
        CssBuilder cellStyle = CssBuilder.EMPTY.paddingLeft(CssSize.em(0.5)).paddingRight(CssSize.em(0.5)).merge(borderStyle);
        CssBuilder coreStyle = CssBuilder.EMPTY.fontWeight("bold");

        Map<String, List<DomContent>> values = Maps.newLinkedHashMap();
        if(aaElements != null && aaElements.variant != null)
        {
            VariantContext variantContext = aaElements.variant;
            List<String> impact = (List<String>) variantContext.getAttribute(IMPACT_KEY);
            String geneName = impact.get(GENE_NAME_IDX);
            String transcriptName = impact.get(TRANSCRIPT_NAME_IDX);
            String aaVariantType = impact.get(AA_VARIANT_TYPE_IDX);
            String hgvs = impact.get(HGVS_INDEX);

            values.put("Gene", List.of(span(geneName)));
            values.put("Transcript", List.of(span(transcriptName)));
            values.put("HGVS", List.of(span(hgvs)));
            values.put("Coding impact", List.of(span(aaVariantType)));
        }

        String filterStr = "PASS";
        if(!filters.isEmpty())
            filterStr = String.join(",", filters);

        List<DomContent> contextElems = Lists.newArrayList();
        contextElems.add(span(mReadContext.leftFlankStr()));
        contextElems.add(span(mReadContext.coreStr()).withStyle(coreStyle.toString()));
        contextElems.add(span(mReadContext.rightFlankStr()));

        String repeatStr = "NO REPEAT";
        if(mReadContext.MaxRepeat != null)
            repeatStr = format("%dx%s", mReadContext.MaxRepeat.Count, mReadContext.MaxRepeat.Bases);

        values.put("Qual", List.of(span(String.valueOf(totalTumorQuality))));
        values.put("Filter", List.of(span(filterStr)));
        values.put("Tier", List.of(span(mVariantTier.name())));
        values.put("CONTEXT", contextElems);
        values.put("NEARBY_INDEL", List.of(span(String.valueOf(nearbyIndel))));
        values.put("REPEAT", List.of(span(repeatStr)));
        values.put("MQF", List.of(span(String.valueOf(mapQualFactor))));

        List<List<TdTag>> rowsElems = Lists.newArrayList();
        for(Map.Entry<String, List<DomContent>> entry : values.entrySet())
        {
            rowsElems.add(List.of(
                    td(entry.getKey()).withStyle(cellStyle.merge(headerStyle).toString()),
                    td().with(entry.getValue()).withStyle(cellStyle.toString())
            ));
        }

        List<DomContent> rows = Lists.newArrayList();
        for(List<TdTag> rowElems : rowsElems)
            rows.add(tr().with(rowElems));

        DomContent table = styledTable(rows, tableStyle);
        return div(table);
    }

    private String getFilename(final String sampleId, @Nullable final String geneName, @Nullable final String variantType)
    {
        String filename = mVariantKey;
        String geneNamePrefix = geneName == null ? "" : geneName + "_";
        String sanitisedVariantType = null;
        if(variantType != null)
        {
            sanitisedVariantType = variantType.replaceAll("&", "_");
            sanitisedVariantType = sanitisedVariantType.replaceAll("\\s", "");
        }

        String variantTypeSuffix = sanitisedVariantType == null ? "" : "_" + sanitisedVariantType;

        if(!filename.startsWith(CHR_PREFIX))
            filename = CHR_PREFIX + filename;

        filename = sampleId + ".sage." + geneNamePrefix + filename;

        SortedSet<String> indexedBasesKeySet = VARIANT_INDEXED_BASES_MAP.get(mVariantKey);
        if(indexedBasesKeySet.size() == 1)
            return filename + variantTypeSuffix + ".html";

        int variantFileNum = indexedBasesKeySet.headSet(mIndexedBasesKey).size() + 1;
        String formatStr = format("%%0%dd", String.valueOf(indexedBasesKeySet.size()).length());
        return filename + "_" + format(formatStr, variantFileNum) + variantTypeSuffix + ".html";
    }

    public void addEvidence(
            final SAMRecord read, @Nullable final FragmentData fragment, final ReadContextMatch matchType,
            @Nullable final QualityScores modifiedQualities)
    {
        ++mReadCount;
        mReadCountByType.put(matchType, mReadCountByType.getOrDefault(matchType, 0) + 1);

        List<ReadEvidenceRecord> records = mReadEvidenceRecordsByType.computeIfAbsent(matchType, k -> Lists.newArrayList());
        if(fragment == null)
        {
            records.add(new ReadEvidenceRecord(read, null, matchType, modifiedQualities, mVariant.Position));
            return;
        }

        BaseRegion firstUnclippedRegion = new BaseRegion(fragment.First.getUnclippedStart(), fragment.First.getUnclippedEnd());
        boolean firstIsVisible = mViewRegion.overlaps(firstUnclippedRegion);

        BaseRegion secondUnclippedRegion = new BaseRegion(fragment.Second.getUnclippedStart(), fragment.Second.getUnclippedEnd());
        boolean secondIsVisible = mViewRegion.overlaps(secondUnclippedRegion);
        if(firstIsVisible && !secondIsVisible)
        {
            records.add(new ReadEvidenceRecord(fragment.First, null, matchType, modifiedQualities, mVariant.Position));
            return;
        }

        if(!firstIsVisible && secondIsVisible)
        {
            records.add(new ReadEvidenceRecord(fragment.Second, null, matchType, modifiedQualities, mVariant.Position));
            return;
        }

        records.add(new ReadEvidenceRecord(read, fragment, matchType, modifiedQualities, mVariant.Position));
    }

    private DomContent renderVariantInfo(final List<String> tumorIds, final List<String> referenceIds)
    {
        CssBuilder strongStyle = CssBuilder.EMPTY.fontWeight("bold");

        StringJoiner infoJoiner = new StringJoiner(" ");
        if(tumorIds.isEmpty())
            infoJoiner.add(referenceIds.get(0));
        else if(referenceIds.isEmpty())
            infoJoiner.add(tumorIds.get(0));
        else
            infoJoiner.add(format("%s[normal: %s]", tumorIds.get(0), referenceIds.get(0)));

        infoJoiner.add(mVariant.chromosome() + ":" + mVariant.Position);
        infoJoiner.add(mVariant.ref() + ">" + mVariant.alt());

        return div(infoJoiner.toString()).withStyle(strongStyle.toString());
    }

    private record AminoAcidElements(String geneRegionLabel, DomContent ref, DomContent alt, @Nullable VariantContext variant) {}

    private record TranscriptNameNode(String name, boolean isVariant, boolean isCanonical) implements Comparable<TranscriptNameNode>
    {

        @Override
        public int compareTo(final TranscriptNameNode o)
        {
            if(isVariant == o.isVariant)
            {
                if(isCanonical == o.isCanonical)
                    return name.compareTo(o.name);

                return isCanonical ? -1 : 1;
            }

            return isVariant ? -1 : 1;
        }
    }

    @Nullable
    private static AminoAcidElements getAminoAcidsElements(final VisConfig config, @Nullable final ReferenceData refData,
            final BaseRegion viewRegion, final SageVariant sageVariant, final RefGenomeSource refGenome)
    {
        if(config.PurpleVcf == null)
        {
            SG_LOGGER.debug("skipping amino acid annotations in vis, since purple vcf not given");
            return null;
        }

        if(refData == null || !refData.geneDataCacheLoaded())
        {
            SG_LOGGER.debug("skipping amino acid annotations in vis, since ensembl data is not loaded");
            return null;
        }

        final SimpleVariant simpleVariant = sageVariant.variant();
        List<VariantContext> vcfVariants;
        try(VcfFileReader vcfFileReader = new VcfFileReader(config.PurpleVcf.toString(), true))
        {
            vcfVariants = vcfFileReader.findVariants(simpleVariant.Chromosome, viewRegion.start(), viewRegion.end());
        }

        if(vcfVariants == null)
            throw new RuntimeException(format("Failed to read purple VCF: %s", config.PurpleVcf));

        if(vcfVariants.isEmpty())
            return null;

        TranscriptNameNode bestTranscriptName = null;
        List<AminoAcidEvent> aaEvents = Lists.newArrayList();
        VariantContext variant = null;
        for(VariantContext variantContext : vcfVariants)
        {
            Object rawImpact = variantContext.getAttribute(IMPACT_KEY);
            if(rawImpact == null)
                continue;

            List<String> impact = (List<String>) rawImpact;
            String transcriptName = impact.get(TRANSCRIPT_NAME_IDX);
            TranscriptAminoAcids transcriptAminoAcids = refData.TransAminoAcidMap.get(transcriptName);
            boolean isCanonical = transcriptAminoAcids.Canonical;

            String hgvs = impact.get(HGVS_INDEX);
            if("p.?".equals(hgvs) || "unknown".equals(hgvs))
            {
                aaEvents = null;
            }
            else if(aaEvents != null && !"".equals(hgvs))
            {
                aaEvents.addAll(AminoAcidEvent.parse(hgvs));
            }

            boolean variantMatches = variantContext.getContig().equals(simpleVariant.chromosome());
            if(variantContext.getStart() != simpleVariant.Position)
                variantMatches = false;

            if(!variantContext.getReference().basesMatch(simpleVariant.Ref))
                variantMatches = false;

            if(!variantContext.getAltAlleleWithHighestAlleleCount().basesMatch(simpleVariant.Alt))
                variantMatches = false;

            TranscriptNameNode newTranscriptName = new TranscriptNameNode(transcriptName, variantMatches, isCanonical);
            if(bestTranscriptName == null || newTranscriptName.compareTo(bestTranscriptName) < 0)
                bestTranscriptName = newTranscriptName;

            if(!variantMatches)
                continue;

            variant = variantContext;
        }

        if(bestTranscriptName == null)
            return null;

        String transcriptName = bestTranscriptName.name;
        TranscriptAminoAcids transcriptAminoAcids = refData.TransAminoAcidMap.get(transcriptName);
        TranscriptData transcriptExons = refData.GeneDataCache.getTranscriptData(transcriptAminoAcids.GeneId, transcriptName);

        if(aaEvents == null)
            aaEvents = Collections.emptyList();

        SvgRender.RenderedGeneData renderedGeneData = renderAminoAcids(
                viewRegion, transcriptExons, transcriptAminoAcids, aaEvents, refGenome, sageVariant);

        String geneRegionLabel = getGeneRegionLabel(transcriptExons, sageVariant.position());
        return new AminoAcidElements(
                transcriptAminoAcids.GeneName + " " + geneRegionLabel,
                rawHtml(renderedGeneData.refSvgCanvas().getSVGElement()),
                rawHtml(renderedGeneData.altSvgCanvas().getSVGElement()),
                variant);
    }

    private List<DomContent> renderReads(boolean isTumor, @Nullable final AminoAcidElements aaElements)
    {
        CssBuilder lightGrayBgStyle = CssBuilder.EMPTY.backgroundColor(Color.LIGHT_GRAY);
        CssBuilder verticalHeaderStyle = CssBuilder.EMPTY.backgroundColor(Color.LIGHT_GRAY).writingMode("vertical-rl");
        CssBuilder headerStyle = CssBuilder.EMPTY.backgroundColor(Color.LIGHT_GRAY).textAlign("right");
        CssBuilder matchTypeBorderStyle = CssBuilder.EMPTY.borderBottom(CssSize.px(2), "solid", Color.BLACK);
        CssBuilder matchTypeStyle = verticalHeaderStyle.merge(matchTypeBorderStyle).textAlign("center").padding(CssSize.px(4));
        CssBuilder tableInfoCellStyle = CssBuilder.EMPTY.fontSizePt((int) Math.round(2.0 / 3.0 * READ_HEIGHT_PX)).fontWeight("bold");
        CssBuilder verticalSpacerDivStyle = CssBuilder.EMPTY.height(CssSize.em(1)).padding(CssSize.ZERO).margin(CssSize.ZERO);

        List<ReadTableColumn> columns = Lists.newArrayList(
                MATE_TYPE_COL, MAP_QUAL_COL, FINAL_QUAL_COL, FINAL_BASE_QUAL_COL, FINAL_MAP_QUAL_COL, SEQ_TECH_BASE_QUAL_COL, ORIENTATION_COL);

        List<DomContent> tableRows = Lists.newArrayList();

        // spacing row
        tableRows.add(tr(td(div().withStyle(verticalSpacerDivStyle.toString())).attr("colspan", columns.size() + 2)));

        // sample row
        String sampleAnnotation = isTumor ? "(tumor)" : "(reference)";
        tableRows.add(tr(td(mSample + " " + sampleAnnotation).attr("colspan", columns.size() + 2)));

        // header row
        List<DomContent> headerCols = Lists.newArrayList();
        headerCols.add(td("Type").withStyle(verticalHeaderStyle.toString()));
        for(ReadTableColumn column : columns)
        {
            headerCols.add(td(column.Header).withStyle(verticalHeaderStyle.toString()));
        }

        headerCols.add(td(rawHtml(renderCoords(READ_HEIGHT_PX, mViewRegion, mVariant.position(), DISPLAY_EVERY_NTH_COORD).getSVGElement())).withStyle(lightGrayBgStyle.toString()));
        DomContent headerRow = tr().with(headerCols);
        tableRows.add(headerRow);

        // amino acids
        if(aaElements != null)
        {
            // ref amino acid row
            DomContent aaRefRow = tr(td("ref").attr("colspan", columns.size() + 1).withStyle(headerStyle.toString()), td(aaElements.ref));
            tableRows.add(aaRefRow);

            // predicted amino acids row
            DomContent aaPredictedRow = tr(
                    td("predicted").attr("colspan", columns.size() + 1).withStyle(headerStyle.toString()), td(aaElements.alt));
            tableRows.add(aaPredictedRow);

            // gene name row
            CssBuilder geneNameStyle = CssBuilder.EMPTY
		    .fontSizePt((int) Math.round(2.0 / 3.0 * READ_HEIGHT_PX)).fontWeight("bold").textAlign("center");
            DomContent geneNameRow = tr(
                    td("gene").attr("colspan", columns.size() + 1).withStyle(headerStyle.toString()),
                    td(aaElements.geneRegionLabel).withStyle(geneNameStyle.toString()));
            tableRows.add(geneNameRow);
        }

        // ref row
        DomContent refRow = tr(td("ref").attr("colspan", columns.size() + 1).withStyle(headerStyle.toString()), td(renderRef()));
        tableRows.add(refRow);

        // context row
        DomContent contextRow =
                tr(td("context").attr("colspan", columns.size() + 1).withStyle(headerStyle.toString()), td(renderContext()));
        tableRows.add(contextRow);

        for(ReadContextMatch matchType : SORTED_MATCH_TYPES)
        {
            List<ReadEvidenceRecord> records = mReadEvidenceRecordsByType.get(matchType);
            if(records == null || records.isEmpty())
            {
                continue;
            }

            Collections.sort(records);
            DomContent typeColContent = rawHtml(
                    format("%s<br>(%d/%d)", matchType.name(), records.size(), mReadCountByType.getOrDefault(matchType, 0)));

            for(int i = 0; i < records.size(); ++i)
            {
                ReadEvidenceRecord record = records.get(i);
                boolean isLastOfType = i == records.size() - 1;

                List<DomContent> cols = Lists.newArrayList();
                if(i == 0)
                {
                    cols.add(td(typeColContent).attr("rowspan", records.size()).withStyle(matchTypeStyle.toString()));
                }

                for(ReadTableColumn column : columns)
                {
                    ContentAndStyle contentAndStyle = column.getContentAndStyle(record);
                    CssBuilder style = tableInfoCellStyle.merge(contentAndStyle.Style);
                    if(isLastOfType)
                    {
                        style = style.merge(matchTypeBorderStyle);
                    }

                    cols.add(contentAndStyle.Content.withStyle(style.toString()));
                }

                if(isLastOfType)
                {
                    cols.add(td(renderRead(record)).withStyle(matchTypeBorderStyle.toString()));
                }
                else
                {
                    cols.add(td(renderRead(record)));
                }

                tableRows.add(tr().with(cols));
            }
        }

        return tableRows;
    }

    private DomContent renderRef()
    {
        return renderBases(mRefViewModel, false, false);
    }

    private DomContent renderContext()
    {
        return renderBases(mContextViewModel, false, true);
    }

    private DomContent renderRead(final ReadEvidenceRecord readEvidence)
    {
        BaseSeqViewModel readViewModel;
        BaseSeqViewModel firstViewModel = null;
        BaseSeqViewModel secondViewModel = null;
        if(readEvidence.Fragment == null)
        {
            readViewModel = BaseSeqViewModel.fromRead(readEvidence.Read);
        }
        else
        {
            firstViewModel = BaseSeqViewModel.fromRead(readEvidence.Fragment.First);
            secondViewModel = BaseSeqViewModel.fromRead(readEvidence.Fragment.Second);
            readViewModel = BaseSeqViewModel.fromConsensusFragment(readEvidence.Read, firstViewModel, secondViewModel);
        }

        SAMRecord firstRead = readEvidence.Read;
        SAMRecord secondRead = null;
        if(readEvidence.Fragment != null)
        {
            firstRead = readEvidence.Fragment.First.getFirstOfPairFlag() ? readEvidence.Fragment.First : readEvidence.Fragment.Second;
            secondRead = readEvidence.Fragment.First.getFirstOfPairFlag() ? readEvidence.Fragment.Second : readEvidence.Fragment.First;
        }

        CssBuilder baseDivStyle = CssBuilder.EMPTY.padding(CssSize.ZERO).margin(CssSize.ZERO);

        DomContent svgEl = renderBases(readViewModel, true, true);
        DomContent svgDiv = div(svgEl).withClass("read-svg").withStyle(baseDivStyle.toString());
        DomContent readInfoDiv = renderReadInfoTable(firstRead, secondRead);

        DomContent containerDiv;
        if(readEvidence.Fragment == null)
        {
            containerDiv = div(svgDiv, readInfoDiv).withStyle(baseDivStyle.toString());
        }
        else
        {
            CssBuilder divStyle = baseDivStyle.display("none");

            DomContent firstSvgEl = renderBases(firstViewModel, true, true);
            DomContent firstSvgDiv = div(firstSvgEl).withClass("read-of-fragment-sgv").withStyle(divStyle.toString());

            DomContent secondSvgEl = renderBases(secondViewModel, true, true);
            DomContent secondSvgDiv = div(secondSvgEl).withClass("read-of-fragment-sgv").withStyle(divStyle.toString());

            containerDiv = div(svgDiv, firstSvgDiv, secondSvgDiv, readInfoDiv).withStyle(baseDivStyle.toString());
        }

        return containerDiv;
    }

    private static SvgRender.RenderedGeneData renderAminoAcids(final BaseRegion viewRegion, final TranscriptData transcriptExons,
            final TranscriptAminoAcids transcriptAminoAcids, final List<AminoAcidEvent> events, final RefGenomeSource refGenome,
            final SageVariant variant)
    {
        boolean posStrand = transcriptExons.posStrand();
        List<GeneRegionViewModel> geneRegions = getGeneRegions(transcriptExons);
        List<GeneRegionViewModel> refViewModels = getRefGeneRegionViewModels(transcriptExons, transcriptAminoAcids, geneRegions);
        List<GeneRegionViewModel> altViewModels = getAltGeneRegionViewModels(
                transcriptExons, transcriptAminoAcids, geneRegions, events, viewRegion, refGenome, variant);
        return renderGeneData(READ_HEIGHT_PX, viewRegion, posStrand, refViewModels, altViewModels);
    }

    private DomContent renderBases(final BaseSeqViewModel bases, boolean shadeQuals, boolean compareToRef)
    {
        SVGGraphics2D svgCanvas =
                renderBaseSeq(READ_HEIGHT_PX, mViewRegion, bases, shadeQuals, mContextBorders, compareToRef ? mRefViewModel : null);
        return rawHtml(svgCanvas.getSVGElement());
    }

    private void downsampleReadEvidenceRecords()
    {
        for(Map.Entry<ReadContextMatch, List<ReadEvidenceRecord>> entry : mReadEvidenceRecordsByType.entrySet())
        {
            ReadContextMatch matchType = entry.getKey();
            List<ReadEvidenceRecord> records = entry.getValue();

            int maxReads = SageVisConstants.MAX_READS_PER_TYPE.get(matchType);

            if(mConfig.MaxSupportReads < 0)
            {
                maxReads = MAX_READ_UPPER_LIMIT;
            }
            else if(mConfig.MaxSupportReads > 0)
            {
                maxReads = min(mConfig.MaxSupportReads, MAX_READ_UPPER_LIMIT);
            }

            List<ReadEvidenceRecord> newRecords = Lists.newArrayList();
            while(!records.isEmpty() && newRecords.size() < maxReads)
            {
                if(records.size() == 1)
                {
                    newRecords.add(records.get(0));
                    records.remove(0);
                    continue;
                }

                int i = ThreadLocalRandom.current().nextInt(records.size());
                if(i < records.size() - 1)
                {
                    ReadEvidenceRecord tmp = records.get(i);
                    records.set(i, records.get(records.size() - 1));
                    records.set(records.size() - 1, tmp);
                }

                newRecords.add(records.get(records.size() - 1));
                records.remove(records.size() - 1);
            }

            mReadEvidenceRecordsByType.put(matchType, newRecords);
        }
    }
}
