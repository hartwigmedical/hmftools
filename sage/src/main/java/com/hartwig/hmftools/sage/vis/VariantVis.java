package com.hartwig.hmftools.sage.vis;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;
import static java.util.Map.entry;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome.CHR_PREFIX;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.CONSENSUS_INFO_DELIM;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.CONSENSUS_READ_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_POSITION;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.UMI_TYPE_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.getMateAlignmentEnd;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.getOrientationString;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.common.NumberEvents.rawNM;
import static com.hartwig.hmftools.sage.vcf.VcfTags.AVG_BASE_QUAL;
import static com.hartwig.hmftools.sage.vcf.VcfTags.AVG_MAP_QUALITY;
import static com.hartwig.hmftools.sage.vcf.VcfTags.FRAG_STRAND_BIAS;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_STRAND_BIAS;
import static com.hartwig.hmftools.sage.vis.ColorUtil.DARK_BLUE;
import static com.hartwig.hmftools.sage.vis.ReadTableColumn.FINAL_QUAL_COL;
import static com.hartwig.hmftools.sage.vis.ReadTableColumn.MAP_QUAL_COL;
import static com.hartwig.hmftools.sage.vis.ReadTableColumn.MATE_TYPE_COL;
import static com.hartwig.hmftools.sage.vis.ReadTableColumn.MOD_BASE_QUAL_COL;
import static com.hartwig.hmftools.sage.vis.ReadTableColumn.MOD_MAP_QUAL_COL;
import static com.hartwig.hmftools.sage.vis.ReadTableColumn.ORIENTATION_COL;
import static com.hartwig.hmftools.sage.vis.ReadTableColumn.RAW_BASE_QUAL_COL;
import static com.hartwig.hmftools.sage.vis.SageVisConstants.BASE_FONT_STYLE;
import static com.hartwig.hmftools.sage.vis.SageVisConstants.DISPLAY_EVERY_NTH_COORD;
import static com.hartwig.hmftools.sage.vis.SageVisConstants.MAX_READ_UPPER_LIMIT;
import static com.hartwig.hmftools.sage.vis.SageVisConstants.READ_HEIGHT_PX;
import static com.hartwig.hmftools.sage.vis.SageVisConstants.VARIANT_INFO_SPACING_SIZE;
import static com.hartwig.hmftools.sage.vis.SvgRender.renderBaseSeq;
import static com.hartwig.hmftools.sage.vis.SvgRender.renderCoords;

import static htsjdk.variant.vcf.VCFConstants.ALLELE_FREQUENCY_KEY;
import static j2html.TagCreator.body;
import static j2html.TagCreator.div;
import static j2html.TagCreator.header;
import static j2html.TagCreator.html;
import static j2html.TagCreator.rawHtml;
import static j2html.TagCreator.script;
import static j2html.TagCreator.span;
import static j2html.TagCreator.table;
import static j2html.TagCreator.td;
import static j2html.TagCreator.tr;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.EnumMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.StringJoiner;
import java.util.TreeSet;
import java.util.concurrent.ThreadLocalRandom;
import java.util.concurrent.atomic.AtomicReference;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.sage.common.ReadContextMatch;
import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.sage.common.RefSequence;
import com.hartwig.hmftools.sage.common.SageVariant;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.common.VariantTier;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;
import com.hartwig.hmftools.sage.quality.QualityScores;
import com.hartwig.hmftools.sage.sync.FragmentData;

import org.jetbrains.annotations.Nullable;
import org.jfree.svg.SVGGraphics2D;

import htsjdk.samtools.SAMRecord;
import j2html.tags.DomContent;
import j2html.tags.specialized.TdTag;

public class VariantVis
{
    private static final AtomicReference<DomContent> JAVASCRIPT = new AtomicReference<>(null);

    private static final DomContent JQUERY_SCRIPT =
            rawHtml("<script src=\"https://ajax.googleapis.com/ajax/libs/jquery/3.7.1/jquery.min.js\"></script>");

    private static final List<ReadContextMatch> SORTED_MATCH_TYPES = Arrays.stream(ReadContextMatch.values())
            .sorted(Comparator.comparingInt(x -> visSortKey(x)))
            .collect(Collectors.toList());

    public static int visSortKey(final ReadContextMatch type)
    {
        switch(type)
        {
            case FULL: return 0;
            case PARTIAL_CORE: return 1;
            case REALIGNED: return 2;
            case CORE: return 3;
            case REF: return 3;
            case NONE:
            default:
                return 5;
        }
    }

    private static Map<String, SortedSet<String>> VARIANT_INDEXED_BASES_MAP = Maps.newConcurrentMap();

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

        int indelSize = mVariant.ref().length() - mVariant.alt().length();
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
        mRefViewModel = BaseSeqViewModel.fromStr(refBases, refPosStart);

        mContextViewModel = BaseSeqViewModel.fromVariant(mReadContext, mVariant.ref(), mVariant.alt());

        StringJoiner indexedBasesKeyBuilder = new StringJoiner("_");
        indexedBasesKeyBuilder.add(String.valueOf(mReadContext.VarIndex));
        indexedBasesKeyBuilder.add(String.valueOf(mReadContext.CoreIndexStart));
        indexedBasesKeyBuilder.add(String.valueOf(mReadContext.CoreIndexEnd));
        indexedBasesKeyBuilder.add(String.valueOf(mReadContext.leftFlankLength()));
        indexedBasesKeyBuilder.add(new String(mReadContext.ReadBases));
        mIndexedBasesKey = indexedBasesKeyBuilder.toString();

        mReadCountByType = Maps.newEnumMap(ReadContextMatch.class);
        mReadCount = 0;

        SortedSet<String> indexedBasesKeySet = VARIANT_INDEXED_BASES_MAP.get(mVariantKey);
        if(indexedBasesKeySet == null)
        {
            indexedBasesKeySet = new TreeSet<>();
            VARIANT_INDEXED_BASES_MAP.put(mVariantKey, indexedBasesKeySet);
        }

        indexedBasesKeySet.add(mIndexedBasesKey);
    }

    private static DomContent styledTable(final List<DomContent> elems, final CssBuilder style)
    {
        return table().with(elems).withStyle(BASE_FONT_STYLE.merge(style).toString());
    }

    public static void writeToHtmlFile(
            final SageVariant sageVariant, final List<String> tumorIds, final List<String> referenceIds, final VisConfig config)
    {
        if(config.PassOnly && !sageVariant.isPassing())
            return;

        List<ReadContextCounter> tumorReadCounters = sageVariant.tumorReadCounters();
        List<ReadContextCounter> refReadCounters = sageVariant.referenceReadCounters();

        // can be null if doesn't meet the configured criteria
        List<VariantVis> tumorVis = tumorReadCounters.stream()
                .map(ReadContextCounter::variantVis).filter(x -> x != null).collect(Collectors.toList());

        List<VariantVis> refVis = refReadCounters.stream()
                .map(ReadContextCounter::variantVis).filter(x -> x != null).collect(Collectors.toList());

        if(tumorVis.isEmpty() && refVis.isEmpty())
            return;

        ReadContextCounter firstCounter = !tumorReadCounters.isEmpty() ? tumorReadCounters.get(0) : refReadCounters.get(0);
        VariantVis firstVis = !tumorVis.isEmpty() ? tumorVis.get(0) : refVis.get(0);

        String sampleId = tumorIds.isEmpty() ? referenceIds.get(0) : tumorIds.get(0);
        String filename = firstVis.getFilename(sampleId);

        int tumorTotalReadCount = tumorVis.stream().mapToInt(x -> x.mReadCount).sum();
        int refTotalReadCount = refVis.stream().mapToInt(x -> x.mReadCount).sum();
        int totalReadCount = tumorTotalReadCount + refTotalReadCount;
        if(totalReadCount == 0)
        {
            SG_LOGGER.info("not writing variant vis file {}, because there are no associated reads", filename);
            return;
        }

        tumorVis.forEach(x -> x.downsampleReadEvidenceRecords());
        refVis.forEach(x -> x.downsampleReadEvidenceRecords());

        Stream<DomContent> tumorReadTableRows = tumorVis.stream().map(x -> x.renderReads(true)).flatMap(x -> x.stream());
        Stream<DomContent> refReadTableRows = refVis.stream().map(x -> x.renderReads(false)).flatMap(x -> x.stream());
        List<DomContent> readTableRows = Stream.concat(tumorReadTableRows, refReadTableRows).collect(Collectors.toList());

        CssBuilder readTableStyle = CssBuilder.EMPTY.borderSpacing(CssSize.ZERO);
        CssBuilder verticalSpacerStyle = CssBuilder.EMPTY.height(CssSize.em(1));

        DomContent readTable = div(styledTable(readTableRows, readTableStyle));
        DomContent verticalSpacer = div().withStyle(verticalSpacerStyle.toString());
        String htmlStr = html(
                header(JQUERY_SCRIPT),
                body(
                        firstVis.renderVariantInfo(
                                sageVariant.totalQuality(),
                                firstCounter.readEdgeDistance().maxAltDistanceFromEdge(), sageVariant.filtersStringSet()),
                        verticalSpacer,
                        renderSampleInfoTable(tumorReadCounters, refReadCounters, tumorIds, referenceIds),
                        readTable,
                        getJavascript()).withStyle(BASE_FONT_STYLE.toString())).render();

        String filePath = Paths.get(firstVis.mConfig.OutputDir, filename).toString();

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

    private static DomContent renderSampleInfoTable(
            final List<ReadContextCounter> tumorReadCounters,
            final List<ReadContextCounter> refReadCounters, final List<String> tumorIds, final List<String> refIds)
    {
        CssBuilder sampleInfoTableStyle = CssBuilder.EMPTY.borderSpacing(CssSize.ZERO);
        CssBuilder firstColCellStyle = CssBuilder.EMPTY.paddingLeft(CssSize.em(0.5)).paddingRight(CssSize.em(0.5));
        CssBuilder cellStyle = firstColCellStyle.textAlign("right");
        CssBuilder firstColHeaderStyle = firstColCellStyle.backgroundColor(Color.LIGHT_GRAY);
        CssBuilder headerStyle = cellStyle.backgroundColor(Color.LIGHT_GRAY);

        List<DomContent> rows = Lists.newArrayList();

        List<String> headers = Lists.newArrayList("SAMPLE", "QUAL", "AD", ALLELE_FREQUENCY_KEY, "DP");
        headers.addAll(SORTED_MATCH_TYPES.stream().map(ReadContextMatch::name).collect(Collectors.toList()));
        headers.addAll(Lists.newArrayList(AVG_BASE_QUAL, AVG_MAP_QUALITY, FRAG_STRAND_BIAS, READ_STRAND_BIAS, "JIT"));

        List<DomContent> headerColumns = Lists.newArrayList();
        for(int i = 0; i < headers.size(); i++)
        {
            CssBuilder style = i == 0 ? firstColHeaderStyle : headerStyle;
            headerColumns.add(td(headers.get(i)).withStyle(style.toString()));
        }

        rows.add(tr().with(headerColumns));

        List<String> allIds = Lists.newArrayList();
        List<ReadContextCounter> allCounters = Lists.newArrayList();
        for (int i = 0; i < tumorReadCounters.size(); i++)
        {
            ReadContextCounter counter = tumorReadCounters.get(i);
            if (counter.variantVis() == null)
                continue;

            allIds.add(tumorIds.get(i));
            allCounters.add(counter);
        }

        for (int i = 0; i < refReadCounters.size(); i++)
        {
            ReadContextCounter counter = refReadCounters.get(i);
            if (counter.variantVis() == null)
                continue;

            allIds.add(refIds.get(i));
            allCounters.add(counter);
        }

        for(int i = 0; i < allCounters.size(); ++i)
        {
            ReadContextCounter counter = allCounters.get(i);
            String sampleId = allIds.get(i);

            int depth = counter.depth();
            int altSupport = counter.altSupport();
            int avgAltMapQuality = altSupport > 0 ? (int) Math.round(counter.altMapQualityTotal() / (double) altSupport) : 0;

            List<TdTag> columnElems = Lists.newArrayList(
                    td(sampleId),
                    td(String.valueOf(counter.tumorQuality())),
                    td(String.valueOf(altSupport)),
                    td(format("%.3f", counter.vaf())),
                    td(String.valueOf(depth)));

            VariantVis variantVis = counter.variantVis();
            for(ReadContextMatch matchType : SORTED_MATCH_TYPES)
            {
                int count = variantVis.mReadCountByType.getOrDefault(matchType, 0);
                columnElems.add(td(String.valueOf(count)));
            }

            columnElems.addAll(Lists.newArrayList(
                    td(String.valueOf((int) counter.averageAltBaseQuality())),
                    td(format("%d", avgAltMapQuality)),
                    td(format("%.2f", counter.fragmentStrandBiasAlt().bias())),
                    td(format("%.2f", counter.readStrandBiasAlt().bias())),
                    td(format("%d-%d", counter.jitter().shortened(), counter.jitter().lengthened()))));

            for(int j = 0; j < columnElems.size(); ++j)
            {
                CssBuilder style = j == 0 ? firstColCellStyle : cellStyle;
                columnElems.set(j, columnElems.get(j).withStyle(style.toString()));
            }

            rows.add(tr().with(columnElems));
        }

        DomContent sampleInfoTable = styledTable(rows, sampleInfoTableStyle);
        return div(sampleInfoTable);
    }

    private static DomContent getJavascript()
    {
        return JAVASCRIPT.updateAndGet((final DomContent currentRef) ->
        {
            if(currentRef != null)
            {
                return currentRef;
            }

            InputStream inputStream = VariantVis.class.getResourceAsStream("/vis/sagevis.js");
            BufferedReader reader = new BufferedReader(new InputStreamReader(inputStream));
            String scriptContent = reader.lines().collect(Collectors.joining("\n"));
            return script(rawHtml(scriptContent)).attr("type", "text/javascript");
        });
    }

    private String getFilename(final String sampleId)
    {
        String filename = mVariantKey;

        if(!filename.startsWith(CHR_PREFIX))
            filename = CHR_PREFIX + filename;

        filename = sampleId + ".sage." + filename;

        SortedSet<String> indexedBasesKeySet = VARIANT_INDEXED_BASES_MAP.get(mVariantKey);
        if(indexedBasesKeySet.size() == 1)
            return filename + ".html";

        int variantFileNum = indexedBasesKeySet.headSet(mIndexedBasesKey).size() + 1;
        String formatStr = format("%%0%dd", String.valueOf(indexedBasesKeySet.size()).length());
        return filename + "_" + format(formatStr, variantFileNum) + ".html";
    }

    public void addEvidence(
            final SAMRecord read, @Nullable final FragmentData fragment, final ReadContextMatch matchType,
            @Nullable final QualityScores modifiedQualities)
    {
        ++mReadCount;
        mReadCountByType.put(matchType, mReadCountByType.getOrDefault(matchType, 0) + 1);

        List<ReadEvidenceRecord> records = mReadEvidenceRecordsByType.get(matchType);
        if(records == null)
        {
            records = Lists.newArrayList();
            mReadEvidenceRecordsByType.put(matchType, records);
        }

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

    private DomContent renderVariantInfo(int totalTumorQuality, int maxDistanceFromEdge, final Set<String> filters)
    {
        CssBuilder horizontalSpacerStyle = CssBuilder.EMPTY.width(VARIANT_INFO_SPACING_SIZE).display("inline-block");
        CssBuilder coreStyle = CssBuilder.EMPTY.fontWeight("bold");

        DomContent horizontalSpacer = div().withStyle(horizontalSpacerStyle.toString());

        String repeatStr = "NO REPEAT";
        if(mReadContext.MaxRepeat != null)
        {
            repeatStr = format("REPEAT = %dx%s", mReadContext.MaxRepeat.Count, mReadContext.MaxRepeat.Bases);
        }

        String filterStr = "FILTER = PASS";
        if(!filters.isEmpty())
            filterStr = "FILTER = " + filters.stream().collect(Collectors.joining(","));

        List<DomContent> contextElems = Lists.newArrayList();
        contextElems.add(span("CONTEXT = "));
        contextElems.add(span(mReadContext.leftFlankStr()));
        contextElems.add(span(mReadContext.coreStr()).withStyle(coreStyle.toString()));
        contextElems.add(span(mReadContext.rightFlankStr()));

        DomContent variantInfoRow = tr(
                td(mVariant.chromosome() + ":" + mVariant.Position),
                td(horizontalSpacer),
                td(mVariant.ref() + " > " + mVariant.alt()),
                td(horizontalSpacer),
                td("TIER = " + mVariantTier.name()),
                td(horizontalSpacer),
                td("QUAL = " + totalTumorQuality),
                td(horizontalSpacer),
                td(repeatStr),
                td(horizontalSpacer),
                td("MED = " + maxDistanceFromEdge),
                td(horizontalSpacer),
                td(filterStr),
                td(horizontalSpacer),
                td().with(contextElems));

        DomContent variantInfoTable = styledTable(Lists.newArrayList(variantInfoRow), CssBuilder.EMPTY);
        return div(variantInfoTable);
    }

    private List<DomContent> renderReads(boolean isTumor)
    {
        CssBuilder lightGrayBgStyle = CssBuilder.EMPTY.backgroundColor(Color.LIGHT_GRAY);
        CssBuilder verticalHeaderStyle = CssBuilder.EMPTY.backgroundColor(Color.LIGHT_GRAY).writingMode("vertical-rl");
        CssBuilder headerStyle = CssBuilder.EMPTY.backgroundColor(Color.LIGHT_GRAY).textAlign("center");
        CssBuilder matchTypeBorderStyle = CssBuilder.EMPTY.borderBottom(CssSize.px(2), "solid", Color.BLACK);
        CssBuilder matchTypeStyle = verticalHeaderStyle.merge(matchTypeBorderStyle).textAlign("center").padding(CssSize.px(4));
        CssBuilder tableInfoCellStyle = CssBuilder.EMPTY.fontSizePt((int) Math.round(2.0 / 3.0 * READ_HEIGHT_PX)).fontWeight("bold");
        CssBuilder verticalSpacerDivStyle = CssBuilder.EMPTY.height(CssSize.em(1)).padding(CssSize.ZERO).margin(CssSize.ZERO);

        List<ReadTableColumn> columns = Lists.newArrayList(
                MATE_TYPE_COL, MAP_QUAL_COL, FINAL_QUAL_COL, MOD_BASE_QUAL_COL, MOD_MAP_QUAL_COL, RAW_BASE_QUAL_COL, ORIENTATION_COL);

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
                    ReadTableColumn.ContentAndStyle contentAndStyle = column.getContentAndStyle(record);
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

    private int getReadNM(final SAMRecord read)
    {
        String chromosome = read.getReferenceName();
        int chromosomeLength = mRefGenCoords.length(chromosome);
        ChrBaseRegion chrRegion = new ChrBaseRegion(chromosome, 1, chromosomeLength);
        RefSequence refSequence = new RefSequence(chrRegion, mRefGenome);
        return rawNM(read, refSequence);
    }

    private DomContent renderReadInfoTable(final SAMRecord firstRead, @Nullable final SAMRecord secondRead)
    {
        CssBuilder baseDivStyle = CssBuilder.EMPTY.padding(CssSize.ZERO).margin(CssSize.ZERO);
        CssBuilder readInfoStyle = baseDivStyle.display("none");

        List<DomContent> readInfoRows = Lists.newArrayList();

        readInfoRows.add(tr(td("Read name:"), td(firstRead.getReadName())));

        String alignmentStr = format("%s:%s-%s", firstRead.getReferenceName(), firstRead.getAlignmentStart(), firstRead.getAlignmentEnd());
        String mateAlignmentStr = "unmapped";
        if(!firstRead.getMateUnmappedFlag())
        {
            String mateChromosome = firstRead.getMateReferenceName();
            int mateAlignmentStart = firstRead.getMateAlignmentStart();
            int mateAlignmentEnd = getMateAlignmentEnd(firstRead);
            String mateAlignmentEndStr = mateAlignmentEnd == NO_POSITION ? "?" : String.valueOf(mateAlignmentEnd);
            mateAlignmentStr = format("%s:%d-%s", mateChromosome, mateAlignmentStart, mateAlignmentEndStr);
        }
        readInfoRows.add(tr(td("Alignment:"), td(alignmentStr + ", " + mateAlignmentStr)));

        String mateCigarStr = firstRead.getStringAttribute(MATE_CIGAR_ATTRIBUTE);
        if(mateCigarStr == null)
        {
            mateCigarStr = "missing";
        }

        readInfoRows.add(tr(td("Cigar:"), td(firstRead.getCigarString() + ", " + mateCigarStr)));

        readInfoRows.add(tr(td("Insert size:"), td(String.valueOf(abs(firstRead.getInferredInsertSize())))));
        readInfoRows.add(tr(td("Orientation:"), td(getOrientationString(firstRead))));

        String firstMapQStr = String.valueOf(firstRead.getMappingQuality());
        String secondMapQStr = secondRead == null ? "" : String.valueOf(secondRead.getMappingQuality());
        String mapQStr = secondRead == null ? firstMapQStr : firstMapQStr + ", " + secondMapQStr;
        readInfoRows.add(tr(td("MapQ:"), td(mapQStr)));

        String firstNumMutationsStr = String.valueOf(getReadNM(firstRead));
        String secondNumMutationsStr = secondRead == null ? "" : String.valueOf(getReadNM(secondRead));
        String numMutationsStr = secondRead == null ? firstNumMutationsStr : firstNumMutationsStr + ", " + secondNumMutationsStr;
        readInfoRows.add(tr(td("NM:"), td(numMutationsStr)));

        String umiTypeStr = firstRead.getStringAttribute(UMI_TYPE_ATTRIBUTE);
        if(umiTypeStr != null)
        {
            readInfoRows.add(tr(td("Dup type:"), td(umiTypeStr)));
        }

        String dupCountStr = "0";
        if(firstRead.hasAttribute(CONSENSUS_READ_ATTRIBUTE))
        {
            dupCountStr = firstRead.getStringAttribute(CONSENSUS_READ_ATTRIBUTE).split(CONSENSUS_INFO_DELIM, 2)[0];
        }

        readInfoRows.add(tr(td("Dup count:"), td(dupCountStr)));

        DomContent readInfoTable = styledTable(readInfoRows, CssBuilder.EMPTY);
        return div(readInfoTable).withClass("read-info").withStyle(readInfoStyle.toString());
    }

    private DomContent renderRead(final ReadEvidenceRecord readEvidence)
    {
        BaseSeqViewModel readViewModel = null;
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