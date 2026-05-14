package com.hartwig.hmftools.sage.vis;

import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.pipeline.MiscToolFiles.SAGE_VIS_PLOT_FILE_EXTENSION;
import static com.hartwig.hmftools.common.pipeline.MiscToolFiles.generateSageVisFilePrefix;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.variant.SageVcfTags.AVG_RECALIBRATED_BASE_QUAL;
import static com.hartwig.hmftools.common.variant.SageVcfTags.MIN_COORDS_COUNT;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_STRAND_BIAS;
import static com.hartwig.hmftools.common.variant.SageVcfTags.UMI_TYPE_COUNTS;
import static com.hartwig.hmftools.common.vis.BaseSeqViewModel.fromStr;
import static com.hartwig.hmftools.common.vis.HtmlUtil.BASE_FONT_STYLE;
import static com.hartwig.hmftools.common.vis.HtmlUtil.getJavascript;
import static com.hartwig.hmftools.common.vis.HtmlUtil.styledTable;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.vcf.VcfTags.AVG_READ_MAP_QUALITY;
import static com.hartwig.hmftools.sage.vcf.VcfTags.AVG_SEQ_TECH_BASE_QUAL;
import static com.hartwig.hmftools.sage.vcf.VcfTags.FRAG_STRAND_BIAS;
import static com.hartwig.hmftools.sage.vis.AminoAcidUtil.getGeneRegionLabel;
import static com.hartwig.hmftools.sage.vis.RenderUtils.renderAminoAcids;
import static com.hartwig.hmftools.sage.vis.RenderUtils.renderReads;
import static com.hartwig.hmftools.sage.vis.RenderUtils.renderVariantInfoTable;
import static com.hartwig.hmftools.sage.vis.VariantVis.VARIANT_INDEXED_BASES_MAP;

import static j2html.TagCreator.body;
import static j2html.TagCreator.div;
import static j2html.TagCreator.html;
import static j2html.TagCreator.rawHtml;
import static j2html.TagCreator.td;
import static j2html.TagCreator.tr;

import java.awt.Color;
import java.io.BufferedWriter;
import java.io.File;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.SortedSet;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.variant.SimpleVariant;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.common.variant.impact.VariantImpactSerialiser;
import com.hartwig.hmftools.common.vis.BaseSeqViewModel;
import com.hartwig.hmftools.common.vis.BaseViewModel;
import com.hartwig.hmftools.common.vis.CssBuilder;
import com.hartwig.hmftools.common.vis.CssSize;
import com.hartwig.hmftools.common.vis.SvgRender;
import com.hartwig.hmftools.sage.ReferenceData;
import com.hartwig.hmftools.sage.common.ReadContextMatch;
import com.hartwig.hmftools.sage.common.SageVariant;
import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;

import org.jetbrains.annotations.Nullable;

import htsjdk.variant.variantcontext.VariantContext;
import j2html.tags.DomContent;
import j2html.tags.specialized.TdTag;

public final class VisFileBuilder
{
    protected static int visSortKey(final ReadContextMatch type)
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

    protected static final List<ReadContextMatch> SORTED_MATCH_TYPES = Arrays.stream(ReadContextMatch.values())
            .sorted(Comparator.comparingInt(VisFileBuilder::visSortKey))
            .toList();

    protected static BaseSeqViewModel contextViewModelFromVariant(final VariantReadContext readContext, final String ref, final String alt)
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

    public static void writeToHtmlFile(
            final SageVariant sageVariant, final List<String> tumorIds, final List<String> referenceIds,
            final VisConfig config, final ReferenceData refData, final RefGenomeSource refGenome)
    {
        VariantVis firstVis = null;

        try
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
            firstVis = !tumorVis.isEmpty() ? tumorVis.get(0) : refVis.get(0);

            AminoAcidElements aaElements = getAminoAcidsElements(config, refData, firstVis.ViewRegion, sageVariant, refGenome);

            String aaVariantType = null;
            String geneName = null;

            if(aaElements != null && aaElements.variant() != null)
            {
                VariantImpact variantImpact = VariantImpactSerialiser.fromVariantContext(aaElements.variant());

                geneName = variantImpact.GeneName;
                aaVariantType = variantImpact.CanonicalEffect;
            }

            String sampleId = tumorIds.isEmpty() ? referenceIds.get(0) : tumorIds.get(0);
            String filename = getFilename(firstVis, sampleId, geneName, aaVariantType);

            int tumorTotalReadCount = tumorVis.stream().mapToInt(x -> x.readCount()).sum();
            int refTotalReadCount = refVis.stream().mapToInt(x -> x.readCount()).sum();
            int totalReadCount = tumorTotalReadCount + refTotalReadCount;
            if(totalReadCount == 0)
            {
                SG_LOGGER.debug("not writing variant vis file {}, because there are no associated reads", filename);
                return;
            }

            tumorVis.forEach(VariantVis::downsampleReadEvidenceRecords);
            refVis.forEach(VariantVis::downsampleReadEvidenceRecords);

            Stream<DomContent> tumorReadTableRows = tumorVis.stream().map(x -> renderReads(x, true, aaElements)).flatMap(Collection::stream);
            Stream<DomContent> refReadTableRows = refVis.stream().map(x -> renderReads(x, false, aaElements)).flatMap(Collection::stream);
            List<DomContent> readTableRows = Stream.concat(tumorReadTableRows, refReadTableRows).collect(Collectors.toList());

            CssBuilder readTableStyle = CssBuilder.EMPTY.borderSpacing(CssSize.ZERO);
            CssBuilder verticalSpacerStyle = CssBuilder.EMPTY.height(CssSize.em(1));

            DomContent readTable = div(styledTable(readTableRows, readTableStyle));

            DomContent verticalSpacer = div().withStyle(verticalSpacerStyle.toString());

            DomContent variantInfo = RenderUtils.renderVariantInfo(firstVis, tumorIds, referenceIds);

            DomContent sampleInfo = styledTable(List.of(tr(
                    td(renderSampleSummaryAndCountsTable(tumorReadCounters, refReadCounters))
                            .withStyle(CssBuilder.EMPTY.verticalAlign("top")
                                    .toString()),
                    td(renderSampleQualAndSiteInfo(tumorReadCounters, refReadCounters, firstCounter))
                            .withStyle(CssBuilder.EMPTY.verticalAlign("top")
                                    .toString()),
                    td(renderVariantInfoTable(
                            firstVis,
                            (int) (-10 * firstCounter.logTqp()),
                            round(round(firstCounter.mapQualFactor() * 10.0d) / 10.0d),
                            sageVariant.nearIndel(),
                            sageVariant.filtersStringSet(),
                            aaElements)).withStyle(CssBuilder.EMPTY.verticalAlign("top").toString()))), CssBuilder.EMPTY);

            String htmlStr = html(
                    body(variantInfo, verticalSpacer, sampleInfo, readTable, getJavascript()).withStyle(BASE_FONT_STYLE.toString())).render();

            String filePath = (new File(firstVis.mConfig.OutputDir, filename)).toString();

            SG_LOGGER.debug("writing variant vis file: {}", filePath);

            BufferedWriter outputWriter = createBufferedWriter(filePath, false);
            outputWriter.write(htmlStr);
            outputWriter.newLine();
            closeBufferedWriter(outputWriter);
        }
        catch(Exception e)
        {
            SG_LOGGER.error("failed to generate variant({}) visualisation:: {}", firstVis.VariantInfo, e.toString());
            System.exit(1);
        }
    }

    protected static String getFilename(
            final VariantVis variant, final String sampleId, @Nullable final String geneName, @Nullable final String variantType)
    {
        SimpleVariant variantInfo = variant.VariantInfo;

        // for reportable variants must conform to
        String variantPrefix = geneName != null ?
                generateSageVisFilePrefix(geneName, variantInfo.chromosome(), variantInfo.position(), variantInfo.Ref, variantInfo.Alt)
                : variant.VariantKey;

        String sanitisedVariantType = null;
        if(variantType != null)
        {
            sanitisedVariantType = variantType.replaceAll("&", "_");
            sanitisedVariantType = sanitisedVariantType.replaceAll("\\s", "");
            sanitisedVariantType = sanitisedVariantType.replaceAll("_variant", "");
        }

        String variantTypeSuffix = sanitisedVariantType == null ? "" : "_" + sanitisedVariantType;

        String filename = sampleId + variantPrefix;

        SortedSet<String> indexedBasesKeySet = VARIANT_INDEXED_BASES_MAP.get(variant.VariantKey);
        if(indexedBasesKeySet.size() == 1)
            return filename + variantTypeSuffix + SAGE_VIS_PLOT_FILE_EXTENSION;

        int variantFileNum = indexedBasesKeySet.headSet(variant.IndexedBasesKey).size() + 1;
        String formatStr = format("%%0%dd", String.valueOf(indexedBasesKeySet.size()).length());
        return filename + "_" + format(formatStr, variantFileNum) + variantTypeSuffix + SAGE_VIS_PLOT_FILE_EXTENSION;
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
            tumorValue = tumorVariantVis == null ? null : String.valueOf(tumorVariantVis.readCountByType().getOrDefault(matchType, 0));
            refValue = refVariantVis == null ? null : String.valueOf(refVariantVis.readCountByType().getOrDefault(matchType, 0));
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

    private static DomContent renderSampleQualAndSiteInfo(
            final List<ReadContextCounter> tumorReadCounters, final List<ReadContextCounter> refReadCounters,
            final ReadContextCounter firstCounter)
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
        values.put("RAW_QUAL", Lists.newArrayList(tumorValue, refValue));

        tumorValue = tumorCounter == null ? null : String.valueOf((int) tumorCounter.averageAltRecalibratedBaseQuality());
        refValue = refCounter == null ? null : String.valueOf((int) refCounter.averageAltRecalibratedBaseQuality());
        values.put(AVG_RECALIBRATED_BASE_QUAL, Lists.newArrayList(tumorValue, refValue));

        tumorValue = tumorCounter == null ? null : String.valueOf((int) tumorCounter.averageAltSeqTechBaseQuality());
        refValue = refCounter == null ? null : String.valueOf((int) refCounter.averageAltSeqTechBaseQuality());
        values.put(AVG_SEQ_TECH_BASE_QUAL, Lists.newArrayList(tumorValue, refValue));

        tumorValue = tumorCounter == null
                ? null
                : String.valueOf(tumorCounter.altSupport() > 0
                        ? (int) round(tumorCounter.altMapQualityTotal() / (double) tumorCounter.altSupport())
                        : 0);
        refValue = refCounter == null
                ? null
                : String.valueOf(
                        refCounter.altSupport() > 0 ? (int) round(refCounter.altMapQualityTotal() / (double) refCounter.altSupport()) : 0);
        values.put(AVG_READ_MAP_QUALITY, Lists.newArrayList(tumorValue, refValue));

        tumorValue = tumorCounter == null ? null : format("%.2f", tumorCounter.fragmentStrandBiasAlt().bias());
        refValue = refCounter == null ? null : format("%.2f", refCounter.fragmentStrandBiasAlt().bias());
        values.put(FRAG_STRAND_BIAS, Lists.newArrayList(tumorValue, refValue));

        tumorValue = tumorCounter == null ? null : format("%.2f", tumorCounter.readStrandBiasAlt().bias());
        refValue = refCounter == null ? null : format("%.2f", refCounter.readStrandBiasAlt().bias());
        values.put(READ_STRAND_BIAS, Lists.newArrayList(tumorValue, refValue));

        tumorValue = tumorCounter == null ? null : format("%d-%d", tumorCounter.jitter().shortened(), tumorCounter.jitter().lengthened());
        refValue = refCounter == null ? null : format("%d-%d", refCounter.jitter().shortened(), refCounter.jitter().lengthened());
        values.put("JIT", Lists.newArrayList(tumorValue, refValue));

        tumorValue = tumorCounter == null ? null : String.valueOf(tumorCounter.fragmentCoords().minCount());
        refValue = refCounter == null ? null : String.valueOf(refCounter.fragmentCoords().minCount());
        values.put(MIN_COORDS_COUNT, Lists.newArrayList(tumorValue, refValue));

        tumorValue = tumorCounter == null ? null : Arrays.toString(tumorCounter.consensusTypeCounts()).replace("[", "").replace("]", "");
        refValue = refCounter == null ? null : Arrays.toString(refCounter.consensusTypeCounts()).replace("[", "").replace("]", "");
        values.put(UMI_TYPE_COUNTS, Lists.newArrayList(tumorValue, refValue));

        tumorValue = tumorCounter == null ? null : format("%.2f", nonAvgEdgeDist);
        refValue = refCounter == null ? null : format("%.2f", altAvgEdgeDist);
        values.put("AED", Lists.newArrayList(tumorValue, refValue));

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

    @Nullable
    private static AminoAcidElements getAminoAcidsElements(
            final VisConfig config, final ReferenceData refData,
            final BaseRegion viewRegion, final SageVariant sageVariant, final RefGenomeSource refGenome)
    {
        if(!config.hasVariantCodingImpacts() || refData == null)
            return null;

        SimpleVariant simpleVariant = sageVariant.variant();

        VariantContext matchedVariantContext = config.getVariantContext(simpleVariant);

        if(matchedVariantContext == null)
            return null;

        List<AminoAcidEvent> aaEvents = Lists.newArrayList();

        VariantImpact variantImpact = VariantImpactSerialiser.fromVariantContext(matchedVariantContext);

        String hgvs = variantImpact.CanonicalHgvsProtein;
        if("p.?".equals(hgvs) || "unknown".equals(hgvs))
        {
            aaEvents = null;
        }
        else if(aaEvents != null && !"".equals(hgvs))
        {
            aaEvents.addAll(AminoAcidEvent.parse(hgvs));
        }

        String transcriptName = variantImpact.CanonicalTranscript;

        TranscriptAminoAcids transcriptAminoAcids = refData.TransAminoAcidMap.get(transcriptName);
        TranscriptData transcriptData = refData.GeneDataCache.getTranscriptData(transcriptAminoAcids.GeneId, transcriptName);

        if(aaEvents == null)
            aaEvents = Collections.emptyList();

        SvgRender.RenderedGeneData renderedGeneData = renderAminoAcids(
                viewRegion, transcriptData, transcriptAminoAcids, aaEvents, refGenome, sageVariant);

        String geneRegionLabel = getGeneRegionLabel(transcriptData, sageVariant.position());

        return new AminoAcidElements(
                transcriptAminoAcids.GeneName + " " + geneRegionLabel,
                rawHtml(renderedGeneData.refSvgCanvas().getSVGElement()),
                rawHtml(renderedGeneData.altSvgCanvas().getSVGElement()),
                matchedVariantContext);
    }
}
