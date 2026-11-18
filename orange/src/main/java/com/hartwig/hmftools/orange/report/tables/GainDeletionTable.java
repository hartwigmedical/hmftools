package com.hartwig.hmftools.orange.report.tables;

import static java.lang.String.format;

import static com.hartwig.hmftools.datamodel.purple.PurpleDriverType.GERMLINE_AMP;
import static com.hartwig.hmftools.datamodel.purple.PurpleDriverType.GERMLINE_DELETION;
import static com.hartwig.hmftools.orange.report.ReportResources.MAX_TABLE_WIDTH;
import static com.hartwig.hmftools.orange.report.ReportResources.NOT_AVAILABLE;
import static com.hartwig.hmftools.orange.report.tables.ExpressionTable.COL_FOLD_CHANGE;
import static com.hartwig.hmftools.orange.report.tables.ExpressionTable.COL_PERCENTILE;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.VALUE_HET;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.VALUE_HOM;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.convertWidths;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatFoldChangeField;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatPercentage;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatSingleDigitDecimal;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatTpmField;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_DRIVER;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_GENE;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_LOCATION;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_RANGE;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_REL_CN;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_TPM;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_TYPE;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.scaleTableWidths;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.createStandardTable;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.createEmptyTable;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.zeroPrefixed;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.purple.PurpleGainDeletion;
import com.hartwig.hmftools.datamodel.purple.PurpleGermlineStatus;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.Cells;

import be.quodlibet.boxable.Cell;
import be.quodlibet.boxable.BaseTable;

import org.apache.pdfbox.pdmodel.PDPage;

import com.hartwig.hmftools.orange.report.DocumentContext;

import java.io.IOException;

import org.apache.logging.log4j.util.Strings;

public final class GainDeletionTable
{
    private static String COL_MIN_CN = "Min CN";
    private static String COL_MAX_CN = "Max CN";
    private static String COL_ARM_CN = "Arm CN";

    public static BaseTable build(final DocumentContext docCtx,
            final String title, float width, final List<PurpleGainDeletion> gainsDels, final ReportResources reportResources,
            boolean hasRna) throws IOException
    {
        if(gainsDels.isEmpty())
        {
            return createEmptyTable(docCtx, title, width, reportResources);
        }

        Cells cells = new Cells(reportResources);

        List<TableHeaderData> headers = Lists.newArrayList();
        headers.add(new TableHeaderData(COL_LOCATION, 0, false));
        headers.add(new TableHeaderData(COL_GENE, 0, false));
        headers.add(new TableHeaderData(COL_TYPE, 0, false));
        headers.add(new TableHeaderData(COL_RANGE, 0, false));
        headers.add(new TableHeaderData(COL_MIN_CN, 0, false));
        headers.add(new TableHeaderData(COL_MAX_CN, 0, false));
        headers.add(new TableHeaderData(COL_REL_CN, 0, false));
        headers.add(new TableHeaderData(COL_ARM_CN, 0, false));

        boolean anyEventsHaveRna = hasRna && gainsDels.stream().anyMatch(x -> x.tpm() != null);

        if(anyEventsHaveRna)
        {
            headers.add(new TableHeaderData(COL_TPM, 0, false));
            headers.add(new TableHeaderData(COL_PERCENTILE, 0, false));
            headers.add(new TableHeaderData(COL_FOLD_CHANGE, 0.6, false));
        }

        headers.add(new TableHeaderData(COL_DRIVER, 0, false));

        List<List<String>> tableRowData = Lists.newArrayList();
        List<Boolean> shadeVariantData = Lists.newArrayList();

        for(PurpleGainDeletion gainDel : sort(gainsDels))
        {
            List<String> rowData = Lists.newArrayList();
            tableRowData.add(rowData);

            rowData.add(gainDel.chromosome() + gainDel.chromosomeBand());
            rowData.add(geneDisplay(gainDel));
            rowData.add(typeDisplay(gainDel));
            rowData.add(rangeDisplay(gainDel));
            rowData.add(formatSingleDigitDecimal(gainDel.minCopyNumber()));
            rowData.add(formatSingleDigitDecimal(gainDel.maxCopyNumber()));
            rowData.add(formatSingleDigitDecimal(gainDel.relativeCopyNumber()));
            rowData.add(formatSingleDigitDecimal(gainDel.armCopyNumber()));

            if(anyEventsHaveRna)
            {
                if(gainDel.tpm() != null)
                {
                    rowData.add(formatTpmField(gainDel.tpm()));
                    rowData.add(formatPercentage(gainDel.tpmPercentile()));
                    rowData.add(formatFoldChangeField(gainDel.tpmFoldChange()));
                }
                else
                {
                    rowData.add(NOT_AVAILABLE);
                    rowData.add(NOT_AVAILABLE);
                    rowData.add(NOT_AVAILABLE);
                }
            }

            rowData.add(gainDel.driver().driverInterpretation().toString());

            boolean isCandidate = gainDel.driver().driverInterpretation() == DriverInterpretation.LOW;
            shadeVariantData.add(isCandidate);
        }

        scaleTableWidths(headers, tableRowData, MAX_TABLE_WIDTH, reportResources);

        int columnCount = headers.size();
        String[] columnTexts = new String[columnCount];
        float[] columnWidths = new float[columnCount];
        float[] columnPercentages = new float[columnCount];

        convertWidths(headers, columnTexts, columnWidths, columnPercentages);

        BaseTable table = createStandardTable(docCtx, title, width, columnWidths, columnTexts, reportResources);

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

    private static List<PurpleGainDeletion> sort(final List<PurpleGainDeletion> gainsAndDels)
    {
        return gainsAndDels.stream().sorted((gainDel1, gainDel2) ->
        {
            double likelihood1 = gainDel1.driver().driverLikelihood();
            double likelihood2 = gainDel2.driver().driverLikelihood();

            int likelihoodCompare = Double.compare(-likelihood1, -likelihood2);

            if(likelihoodCompare != 0)
            {
                return likelihoodCompare;
            }

            String location1 = zeroPrefixed(gainDel1.chromosome() + gainDel1.chromosomeBand());
            String location2 = zeroPrefixed(gainDel2.chromosome() + gainDel2.chromosomeBand());

            if(location1.equals(location2))
            {
                return gainDel1.gene().compareTo(gainDel2.gene());
            }
            else
            {
                return location1.compareTo(location2);
            }
        }).collect(Collectors.toList());
    }

    private static String rangeDisplay(final PurpleGainDeletion gainDel)
    {
        if(gainDel.exonStart() == null && gainDel.exonEnd() == null)
        {
            return gainDel.geneRange();
        }

        return format("Exon %d - Exon %d", gainDel.exonStart(), gainDel.exonEnd());
    }

    private static String typeDisplay(final PurpleGainDeletion gainDel)
    {
        if(gainDel.driver().type() == GERMLINE_AMP)
        {
            return "AMP";
        }
        else if(gainDel.driver().type() == GERMLINE_DELETION)
        {
            return gainDel.germlineAmpDelFields().germlineStatus() == PurpleGermlineStatus.HOM_DELETION ? VALUE_HOM : VALUE_HET;
        }

        return gainDel.driver().type().toString();
    }

    private static String geneDisplay(final PurpleGainDeletion gainDel)
    {
        String addon = Strings.EMPTY;
        if(!gainDel.isCanonical())
        {
            addon = " (alt)";
        }
        return gainDel.gene() + addon;
    }
}
