package com.hartwig.hmftools.orange.report.tables;

import static java.lang.String.format;

import static com.hartwig.hmftools.orange.algo.OrangeConstants.isCandidateLikelihood;
import static com.hartwig.hmftools.orange.report.ReportResources.MAX_TABLE_WIDTH;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_COPIES;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_GENE;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_HGVS;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_POSITION;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.VALUE_NO;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.VALUE_YES;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.convertWidths;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatPercentageField;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatTwoDigitDecimal;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_AF;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_BIALLELIC;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_CL;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_DRIVER;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_HOTSPOT;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_SL;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_DP;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_RNA;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.scaleTableWidths;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.createStandardTable;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.createEmptyTable;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.datamodel.common.AllelicDepth;
import com.hartwig.hmftools.datamodel.purple.PurpleTranscriptImpact;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.Cells;

import be.quodlibet.boxable.Cell;
import be.quodlibet.boxable.BaseTable;

import org.apache.pdfbox.pdmodel.PDPage;

import com.hartwig.hmftools.orange.report.DocumentContext;

import java.io.IOException;

public final class SomaticVariantTable
{
    public static BaseTable build(final DocumentContext docCtx,
            final String title, float width, final List<PurpleVariant> variants, final ReportResources reportResources,
            boolean tumorOnly, boolean hasRna) throws IOException
    {
        if(variants.isEmpty())
        {
            return createEmptyTable(docCtx, title, width, reportResources);
        }

        List<TableHeaderData> headers = Lists.newArrayList();

        headers.add(new TableHeaderData(COL_GENE, 0, false));
        headers.add(new TableHeaderData(COL_POSITION, 0, false));
        headers.add(new TableHeaderData(COL_HGVS, 0, true));
        headers.add(new TableHeaderData(COL_AF, 0, false));
        headers.add(new TableHeaderData(COL_DP, 0, false));
        headers.add(new TableHeaderData(COL_COPIES, 0, false));
        headers.add(new TableHeaderData(COL_HOTSPOT, 0, false));
        headers.add(new TableHeaderData(COL_BIALLELIC, 0, false));
        headers.add(new TableHeaderData(COL_CL, 0, false));

        if(tumorOnly)
            headers.add(new TableHeaderData(COL_SL, 0, false));

        if(hasRna)
            headers.add(new TableHeaderData(COL_RNA, 0, false));

        headers.add(new TableHeaderData(COL_DRIVER, 0, false));

        List<List<String>> tableRowData = Lists.newArrayList();
        List<Boolean> shadeVariantData = Lists.newArrayList();

        for(PurpleVariant variant : sort(variants))
        {
            List<PurpleTranscriptImpact> transcriptImpacts = Lists.newArrayList(variant.canonicalImpact());
            transcriptImpacts.addAll(variant.otherImpacts());

            for(PurpleTranscriptImpact transcriptImpact : transcriptImpacts)
            {
                List<String> rowData = Lists.newArrayList();
                tableRowData.add(rowData);

                rowData.add(variant.gene());
                rowData.add(locationDisplay(variant));
                rowData.add(hgvsDisplay(transcriptImpact));
                rowData.add(formatTwoDigitDecimal(variant.tumorDepth().alleleFrequency()));
                rowData.add(String.valueOf(variant.tumorDepth().totalReadCount()));
                rowData.add(copyNumberDisplay(variant));
                rowData.add(hotspotDisplay(variant));
                rowData.add(formatPercentageField(variant.biallelicProbability()));

                rowData.add(clonalLikelihood(variant));

                if(tumorOnly)
                {
                    rowData.add(variant.somaticLikelihood().toString());
                }

                if(hasRna)
                {
                    rowData.add(rnaDisplay(variant));
                }

                rowData.add(formatPercentageField(variant.driverLikelihood()));

                shadeVariantData.add(isCandidateLikelihood(variant.driverLikelihood()));
            }
        }

        scaleTableWidths(headers, tableRowData, MAX_TABLE_WIDTH, reportResources);

        int columnCount = headers.size();
        String[] columnTexts = new String[columnCount];
        float[] columnWidths = new float[columnCount];
        float[] columnPercentages = new float[columnCount];

        convertWidths(headers, columnTexts, columnWidths, columnPercentages);

        BaseTable table = createStandardTable(docCtx, title, width, columnWidths, columnTexts, reportResources);

        Cells cells = new Cells(reportResources);

        for(int r = 0; r < tableRowData.size(); ++r)
        {
            List<String> rowData = tableRowData.get(r);
            List<Cell<PDPage>> createdCells = cells.addRow(table, columnPercentages, rowData);

            if(shadeVariantData.get(r))
            {
                reportResources.shadeCandidateCells(createdCells);
            }
        }

        return table;
    }

    protected static List<PurpleVariant> sort(final List<PurpleVariant> variants)
    {
        return variants.stream().sorted((variant1, variant2) ->
        {
            int driverCompare = Double.compare(variant2.driverLikelihood(), variant1.driverLikelihood());
            if(driverCompare != 0)
            {
                return driverCompare;
            }

            int geneCompare = variant1.gene().compareTo(variant2.gene());
            if(geneCompare != 0)
            {
                return geneCompare;
            }

            return 0;
        }).collect(Collectors.toList());
    }

    protected static String locationDisplay(final PurpleVariant variant)
    {
        boolean phased = variant.localPhaseSets() != null && !variant.localPhaseSets().isEmpty();
        String locationStr = format("%s:%d", variant.chromosome(), variant.position());
        return phased ? format("%s *", locationStr) : locationStr;
    }

    protected static String hgvsDisplay(final PurpleTranscriptImpact transcriptImpact)
    {
        if(transcriptImpact.hgvsProteinImpact().isEmpty())
        {
            return transcriptImpact.hgvsCodingImpact();
        }

        return format("%s [%s]", transcriptImpact.hgvsCodingImpact(), transcriptImpact.hgvsProteinImpact());
    }

    protected static String copyNumberDisplay(final PurpleVariant variant)
    {
        return format("%.1f of %.1f", variant.variantCopyNumber(), variant.adjustedCopyNumber());
    }

    protected static String clonalLikelihood(final PurpleVariant variant)
    {
        return formatPercentageField(1 - variant.subclonalLikelihood());
    }

    protected static String hotspotDisplay(final PurpleVariant variant)
    {
        switch(variant.hotspot())
        {
            case HOTSPOT:
                return VALUE_YES;
            case NEAR_HOTSPOT:
                return "Near";
            default:
                return VALUE_NO;
        }
    }

    protected static String rnaDisplay(final PurpleVariant variant)
    {
        AllelicDepth rnaDepth = variant.rnaDepth();

        if(rnaDepth == null)
        {
            return ReportResources.NOT_AVAILABLE;
        }

        return format("%d/%d", rnaDepth.alleleReadCount(), rnaDepth.totalReadCount());
    }
}

