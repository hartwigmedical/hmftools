package com.hartwig.hmftools.orange.report.tables;

import static java.lang.String.format;

import static com.hartwig.hmftools.orange.algo.OrangeConstants.isCandidateLikelihood;
import static com.hartwig.hmftools.orange.report.ReportResources.MAX_TABLE_WIDTH;
import static com.hartwig.hmftools.orange.report.tables.SomaticVariantTable.copyNumberDisplay;
import static com.hartwig.hmftools.orange.report.tables.SomaticVariantTable.hgvsDisplay;
import static com.hartwig.hmftools.orange.report.tables.SomaticVariantTable.hotspotDisplay;
import static com.hartwig.hmftools.orange.report.tables.SomaticVariantTable.locationDisplay;
import static com.hartwig.hmftools.orange.report.tables.SomaticVariantTable.sort;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_COPIES;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_GENE;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_HGVS;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.VALUE_NO;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.VALUE_YES;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.convertWidths;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatTwoDigitDecimal;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_AF;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_BIALLELIC;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_HOTSPOT;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_VARIANT;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_DP;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_DRIVER;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_ZYGOSITY;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.VALUE_HET;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.VALUE_HOM;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.scaleTableWidths;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.createStandardTable;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.createEmptyTable;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.purple.PurpleGenotypeStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleTranscriptImpact;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.Cells;

import be.quodlibet.boxable.Cell;
import be.quodlibet.boxable.BaseTable;

import org.apache.pdfbox.pdmodel.PDPage;

import com.hartwig.hmftools.orange.report.DocumentContext;

import java.io.IOException;

public final class GermlineVariantTable
{
    public static final String COL_CLINVAR = "Clinvar";
    public static final String COL_GNOMAD = "Gnomad";

    public static BaseTable build(final DocumentContext docCtx, final String title, float width, final List<PurpleVariant> variants,
            final ReportResources reportResources) throws IOException
    {
        if(variants.isEmpty())
        {
            return createEmptyTable(docCtx, title, width, reportResources);
        }

        List<TableHeaderData> headers = Lists.newArrayList();

        headers.add(new TableHeaderData(COL_GENE, 0, false));
        headers.add(new TableHeaderData(COL_VARIANT, 0, false));
        headers.add(new TableHeaderData(COL_HGVS, 0, true));
        headers.add(new TableHeaderData(COL_ZYGOSITY, 0, false));
        headers.add(new TableHeaderData(COL_AF, 0, false));
        headers.add(new TableHeaderData(COL_DP, 0, false));
        headers.add(new TableHeaderData(COL_COPIES, 0, false));
        headers.add(new TableHeaderData(COL_HOTSPOT, 0, false));
        headers.add(new TableHeaderData(COL_BIALLELIC, 0, false));
        headers.add(new TableHeaderData(COL_GNOMAD, 0, false));
        headers.add(new TableHeaderData(COL_CLINVAR, 0, false));
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
                rowData.add(simplifiedDisplay(variant.genotypeStatus()));
                rowData.add(formatTwoDigitDecimal(variant.tumorDepth().alleleFrequency()));
                rowData.add(String.valueOf(variant.tumorDepth().totalReadCount()));
                rowData.add(copyNumberDisplay(variant));
                rowData.add(hotspotDisplay(variant));
                rowData.add(variant.biallelic() ? VALUE_YES : VALUE_NO);
                rowData.add(format("%4.2e", variant.gnomadFrequency()));
                rowData.add(variant.clinvarPathogenicity());

                DriverInterpretation driverInterpretation = DriverInterpretation.interpret(variant.driverLikelihood());
                rowData.add(driverInterpretation.toString());

                boolean isCandidate = isCandidateLikelihood(variant.driverLikelihood());
                shadeVariantData.add(isCandidate);
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

    private static String simplifiedDisplay(PurpleGenotypeStatus genotypeStatus)
    {
        return switch(genotypeStatus)
        {
            case HOM_REF, HOM_ALT -> VALUE_HOM;
            case HET -> VALUE_HET;
            case UNKNOWN -> "UNKNOWN";
        };
    }
}
