package com.hartwig.hmftools.sage.vis;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.vis.HtmlUtil.renderReadInfoTable;
import static com.hartwig.hmftools.common.vis.HtmlUtil.styledTable;
import static com.hartwig.hmftools.common.vis.SvgRender.renderBaseSeq;
import static com.hartwig.hmftools.common.vis.SvgRender.renderCoords;
import static com.hartwig.hmftools.common.vis.SvgRender.renderGeneData;
import static com.hartwig.hmftools.sage.vis.AminoAcidUtil.getAltGeneRegionViewModels;
import static com.hartwig.hmftools.sage.vis.AminoAcidUtil.getGeneRegions;
import static com.hartwig.hmftools.sage.vis.AminoAcidUtil.getRefGeneRegionViewModels;
import static com.hartwig.hmftools.sage.vis.ReadTableColumn.FINAL_BASE_QUAL_COL;
import static com.hartwig.hmftools.sage.vis.ReadTableColumn.FINAL_MAP_QUAL_COL;
import static com.hartwig.hmftools.sage.vis.ReadTableColumn.FINAL_QUAL_COL;
import static com.hartwig.hmftools.sage.vis.ReadTableColumn.MAP_QUAL_COL;
import static com.hartwig.hmftools.sage.vis.ReadTableColumn.MATE_TYPE_COL;
import static com.hartwig.hmftools.sage.vis.ReadTableColumn.ORIENTATION_COL;
import static com.hartwig.hmftools.sage.vis.ReadTableColumn.SEQ_TECH_BASE_QUAL_COL;
import static com.hartwig.hmftools.sage.vis.SageVisConstants.DISPLAY_EVERY_NTH_COORD;
import static com.hartwig.hmftools.sage.vis.SageVisConstants.READ_HEIGHT_PX;
import static com.hartwig.hmftools.sage.vis.VisFileBuilder.SORTED_MATCH_TYPES;

import static j2html.TagCreator.div;
import static j2html.TagCreator.rawHtml;
import static j2html.TagCreator.span;
import static j2html.TagCreator.td;
import static j2html.TagCreator.tr;

import java.awt.Color;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.variant.CommonVcfTags;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.common.variant.impact.VariantImpactSerialiser;
import com.hartwig.hmftools.common.vis.BaseSeqViewModel;
import com.hartwig.hmftools.common.vis.CssBuilder;
import com.hartwig.hmftools.common.vis.CssSize;
import com.hartwig.hmftools.common.vis.GeneRegionViewModel;
import com.hartwig.hmftools.common.vis.SvgRender;
import com.hartwig.hmftools.sage.common.ReadContextMatch;
import com.hartwig.hmftools.sage.common.SageVariant;

import org.jetbrains.annotations.Nullable;
import org.jfree.svg.SVGGraphics2D;

import htsjdk.samtools.SAMRecord;
import j2html.tags.DomContent;
import j2html.tags.specialized.TdTag;

public class RenderUtils
{
    protected static DomContent renderVariantInfoTable(
            final VariantVis variant, int totalTumorQuality, double mapQualFactor, boolean nearbyIndel, final Set<String> filters,
            @Nullable final AminoAcidElements aaElements)
    {
        CssBuilder borderStyle = CssBuilder.EMPTY.border(CssSize.px(1.0), "solid", Color.BLACK).borderCollapse("collapse");
        CssBuilder tableStyle = borderStyle;
        CssBuilder headerStyle = CssBuilder.EMPTY.fontWeight("bold").textAlign("center").backgroundColor(Color.LIGHT_GRAY);
        CssBuilder cellStyle = CssBuilder.EMPTY.paddingLeft(CssSize.em(0.5)).paddingRight(CssSize.em(0.5)).merge(borderStyle);
        CssBuilder coreStyle = CssBuilder.EMPTY.fontWeight("bold");

        Map<String, List<DomContent>> values = Maps.newLinkedHashMap();

        if(aaElements != null && aaElements.variant() != null)
        {
            VariantImpact variantImpact = VariantImpactSerialiser.fromVariantContext(aaElements.variant());

            values.put("Gene", List.of(span(variantImpact.GeneName)));
            values.put("Transcript", List.of(span(variantImpact.CanonicalTranscript)));
            values.put("HGVS", List.of(span(variantImpact.CanonicalHgvsProtein)));
            values.put("Coding impact", List.of(span(variantImpact.CanonicalEffect)));
        }

        String filterStr = CommonVcfTags.PASS_FILTER;
        if(!filters.isEmpty())
            filterStr = String.join(",", filters);

        List<DomContent> contextElems = Lists.newArrayList();
        contextElems.add(span(variant.ReadContext.leftFlankStr()));
        contextElems.add(span(variant.ReadContext.coreStr()).withStyle(coreStyle.toString()));
        contextElems.add(span(variant.ReadContext.rightFlankStr()));

        String repeatStr = "NO REPEAT";
        if(variant.ReadContext.MaxRepeat != null)
            repeatStr = format("%dx%s", variant.ReadContext.MaxRepeat.Count, variant.ReadContext.MaxRepeat.Bases);

        values.put("Qual", List.of(span(String.valueOf(totalTumorQuality))));
        values.put("Filter", List.of(span(filterStr)));
        values.put("Tier", List.of(span(variant.Tier.name())));
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

    public static List<DomContent> renderReads(final VariantVis variant, boolean isTumor, @Nullable final AminoAcidElements aaElements)
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
        tableRows.add(tr(td(variant.SampleId + " " + sampleAnnotation).attr("colspan", columns.size() + 2)));

        // header row
        List<DomContent> headerCols = Lists.newArrayList();
        headerCols.add(td("Type").withStyle(verticalHeaderStyle.toString()));
        for(ReadTableColumn column : columns)
        {
            headerCols.add(td(column.Header).withStyle(verticalHeaderStyle.toString()));
        }

        headerCols.add(td(rawHtml(renderCoords(
                READ_HEIGHT_PX, variant.ViewRegion, variant.VariantInfo.position(),
                DISPLAY_EVERY_NTH_COORD).getSVGElement())).withStyle(lightGrayBgStyle.toString()));
        DomContent headerRow = tr().with(headerCols);
        tableRows.add(headerRow);

        // amino acids
        if(aaElements != null)
        {
            // ref amino acid row
            DomContent aaRefRow = tr(
                    td("ref").attr("colspan", columns.size() + 1).withStyle(headerStyle.toString()), td(aaElements.ref()));
            tableRows.add(aaRefRow);

            // predicted amino acids row
            DomContent aaPredictedRow = tr(
                    td("predicted").attr("colspan", columns.size() + 1).withStyle(headerStyle.toString()),
                    td(aaElements.alt()));
            tableRows.add(aaPredictedRow);

            // gene name row
            CssBuilder geneNameStyle = CssBuilder.EMPTY
                    .fontSizePt((int) Math.round(2.0 / 3.0 * READ_HEIGHT_PX)).fontWeight("bold").textAlign("center");
            DomContent geneNameRow = tr(
                    td("gene").attr("colspan", columns.size() + 1).withStyle(headerStyle.toString()),
                    td(aaElements.geneRegionLabel()).withStyle(geneNameStyle.toString()));
            tableRows.add(geneNameRow);
        }

        // ref row
        DomContent refRow = tr(
                td("ref").attr("colspan", columns.size() + 1).withStyle(headerStyle.toString()), td(renderRef(variant)));
        tableRows.add(refRow);

        // context row
        DomContent contextRow = tr(
                td("context").attr("colspan", columns.size() + 1).withStyle(headerStyle.toString()),
                td(renderContext(variant)));
        tableRows.add(contextRow);

        for(ReadContextMatch matchType : SORTED_MATCH_TYPES)
        {
            List<ReadEvidenceRecord> records = variant.readEvidenceRecordsByType().get(matchType);
            if(records == null || records.isEmpty())
            {
                continue;
            }

            Collections.sort(records);
            DomContent typeColContent = rawHtml(
                    format("%s<br>(%d/%d)", matchType.name(), records.size(), variant.readCountByType().getOrDefault(matchType, 0)));

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
                    cols.add(td(renderRead(variant, record)).withStyle(matchTypeBorderStyle.toString()));
                }
                else
                {
                    cols.add(td(renderRead(variant, record)));
                }

                tableRows.add(tr().with(cols));
            }
        }

        return tableRows;
    }

    private static DomContent renderRef(final VariantVis variant)
    {
        return renderBases(variant, variant.RefViewModel, false, false);
    }

    private static DomContent renderContext(final VariantVis variant)
    {
        return renderBases(variant, variant.ContextViewModel, false, true);
    }

    private static DomContent renderRead(final VariantVis variant, final ReadEvidenceRecord readEvidence)
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

        DomContent svgEl = renderBases(variant, readViewModel, true, true);
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

            DomContent firstSvgEl = renderBases(variant, firstViewModel, true, true);
            DomContent firstSvgDiv = div(firstSvgEl).withClass("read-of-fragment-sgv").withStyle(divStyle.toString());

            DomContent secondSvgEl = renderBases(variant, secondViewModel, true, true);
            DomContent secondSvgDiv = div(secondSvgEl).withClass("read-of-fragment-sgv").withStyle(divStyle.toString());

            containerDiv = div(svgDiv, firstSvgDiv, secondSvgDiv, readInfoDiv).withStyle(baseDivStyle.toString());
        }

        return containerDiv;
    }

    protected static SvgRender.RenderedGeneData renderAminoAcids(
            final BaseRegion viewRegion, final TranscriptData transcriptExons,
            final TranscriptAminoAcids transcriptAminoAcids, final List<AminoAcidEvent> events, final RefGenomeSource refGenome,
            final SageVariant sageVariant)
    {
        boolean posStrand = transcriptExons.posStrand();

        List<GeneRegionViewModel> geneRegions = getGeneRegions(transcriptExons);

        List<GeneRegionViewModel> refViewModels = getRefGeneRegionViewModels(transcriptExons, transcriptAminoAcids, geneRegions);

        List<GeneRegionViewModel> altViewModels = getAltGeneRegionViewModels(
                transcriptExons, transcriptAminoAcids, geneRegions, events, viewRegion, refGenome, sageVariant);

        return renderGeneData(READ_HEIGHT_PX, viewRegion, posStrand, refViewModels, altViewModels);
    }

    private static DomContent renderBases(final VariantVis variant, final BaseSeqViewModel bases, boolean shadeQuals, boolean compareToRef)
    {
        SVGGraphics2D svgCanvas =
                renderBaseSeq(READ_HEIGHT_PX, variant.ViewRegion, bases, shadeQuals, variant.ContextBorders,
                        compareToRef ? variant.RefViewModel : null);

        return rawHtml(svgCanvas.getSVGElement());
    }

    protected static DomContent renderVariantInfo(final VariantVis variant, final List<String> tumorIds, final List<String> referenceIds)
    {
        CssBuilder strongStyle = CssBuilder.EMPTY.fontWeight("bold");

        StringJoiner infoJoiner = new StringJoiner(" ");
        if(tumorIds.isEmpty())
            infoJoiner.add(referenceIds.get(0));
        else if(referenceIds.isEmpty())
            infoJoiner.add(tumorIds.get(0));
        else
            infoJoiner.add(format("%s[normal: %s]", tumorIds.get(0), referenceIds.get(0)));

        infoJoiner.add(variant.VariantInfo.chromosome() + ":" + variant.VariantInfo.Position);
        infoJoiner.add(variant.VariantInfo.ref() + ">" + variant.VariantInfo.alt());

        return div(infoJoiner.toString()).withStyle(strongStyle.toString());
    }
}
