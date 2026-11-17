package com.hartwig.hmftools.orange.report.tables;

import static java.lang.String.format;

import static com.hartwig.hmftools.orange.report.ReportResources.MAX_TABLE_WIDTH;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_JUNCTIONS;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_RNA_FRAGS;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_TYPE;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.convertWidths;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatSingleDigitDecimal;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_DRIVER;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_FUSION;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_JCN;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.scaleTableWidths;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.createStandardTable;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.createEmptyTable;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.linx.FusionPhasedType;
import com.hartwig.hmftools.datamodel.linx.LinxFusion;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.Cells;

import be.quodlibet.boxable.BaseTable;

import com.hartwig.hmftools.orange.report.DocumentContext;

import java.io.IOException;

import com.google.common.collect.Lists;

public final class DnaFusionTable
{
    private static final String COL_PHASING = "Phasing";

    public static BaseTable build(final DocumentContext docCtx,
            final String title, float width, final List<LinxFusion> fusions, final ReportResources reportResources) throws IOException
    {
        if(fusions.isEmpty())
        {
            return createEmptyTable(docCtx, title, width, reportResources);
        }

        List<TableHeaderData> headers = Lists.newArrayList();

        headers.add(new TableHeaderData(COL_FUSION, 0, false));
        headers.add(new TableHeaderData(COL_JUNCTIONS, 0, false));
        headers.add(new TableHeaderData(COL_JCN, 0, false));
        headers.add(new TableHeaderData(COL_PHASING, 0, false));
        headers.add(new TableHeaderData(COL_TYPE, 0, false));

        boolean hasRna = fusions.stream().anyMatch(x -> x.rnaSupport() != null);

        if(hasRna)
            headers.add(new TableHeaderData(COL_RNA_FRAGS, 0, false));

        headers.add(new TableHeaderData(COL_DRIVER, 0, false));

        List<List<String>> tableRowData = Lists.newArrayList();

        for(LinxFusion fusion : sortLinxFusions(fusions))
        {
            List<String> rowData = Lists.newArrayList();
            tableRowData.add(rowData);

            rowData.add(fusionDisplay(fusion));
            rowData.add(transcriptJunctions(fusion));
            rowData.add(formatSingleDigitDecimal(fusion.junctionCopyNumber()));
            rowData.add(display(fusion.phased()));
            rowData.add(fusion.reportedType().toString());

            if(hasRna)
            {
                rowData.add(String.valueOf(fusion.rnaSupport().alleleReadCount()));
            }

            rowData.add(fusion.driverInterpretation().toString());
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
            cells.addRow(table, columnPercentages, rowData);
        }

        return table;
    }

    private static String fusionDisplay(final LinxFusion fusion)
    {
        return format("%s::%s", fusion.geneUp(), fusion.geneDown());
    }

    private static String transcriptJunctions(final LinxFusion fusion)
    {
        return format("%s (%s) - %s (%s)",
                fusion.contextUp(), fusion.transcriptUp(), fusion.contextDown(), fusion.transcriptDown());
    }

    private static String display(FusionPhasedType fusionPhasedType)
    {
        return switch(fusionPhasedType)
        {
            case INFRAME -> "Inframe";
            case SKIPPED_EXONS -> "Skipped exons";
            case OUT_OF_FRAME -> "Out of frame";
        };
    }

    private static List<LinxFusion> sortLinxFusions(final List<LinxFusion> fusions)
    {
        return fusions.stream().sorted((fusion1, fusion2) ->
        {
            if(fusion1.driverInterpretation() == fusion2.driverInterpretation())
            {
                if(fusion1.geneUp().equals(fusion2.geneUp()))
                {
                    return fusion1.geneDown().compareTo(fusion2.geneDown());
                }
                else
                {
                    return fusion1.geneUp().compareTo(fusion2.geneUp());
                }
            }
            else
            {
                return fusion1.driverInterpretation() == DriverInterpretation.HIGH ? -1 : 1;
            }
        }).collect(Collectors.toList());
    }
}
