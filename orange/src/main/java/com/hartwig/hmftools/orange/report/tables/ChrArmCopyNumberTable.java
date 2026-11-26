package com.hartwig.hmftools.orange.report.tables;

import static com.hartwig.hmftools.orange.report.ReportResources.MAX_TABLE_WIDTH;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_DRIVER;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.convertWidths;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatSingleDigitDecimal;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_CHR;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_CN;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_REL_CN;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_TYPE;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.scaleTableWidths;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.createStandardTable;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.createEmptyTable;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.purple.PurpleChrArmCopyNumber;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.Cells;

import be.quodlibet.boxable.Cell;
import be.quodlibet.boxable.BaseTable;

import org.apache.pdfbox.pdmodel.PDPage;

import com.hartwig.hmftools.orange.report.DocumentContext;

import java.io.IOException;

public final class ChrArmCopyNumberTable
{
    private static final String COL_ARM = "Arm";

    public static BaseTable build(final DocumentContext docCtx,
            final String title, float width, final List<PurpleChrArmCopyNumber> chrArmCopyNumbers, final ReportResources reportResources)
            throws IOException
    {
        if(chrArmCopyNumbers.isEmpty())
        {
            return createEmptyTable(docCtx, title, width, reportResources);
        }

        List<TableHeaderData> headers = Lists.newArrayList();

        headers.add(new TableHeaderData(COL_CHR, 0, false));
        headers.add(new TableHeaderData(COL_ARM, 0, false));
        headers.add(new TableHeaderData(COL_TYPE, 0, false));
        headers.add(new TableHeaderData(COL_CN, 0, false));
        headers.add(new TableHeaderData(COL_REL_CN, 0, false));
        headers.add(new TableHeaderData(COL_DRIVER, 0, false));

        List<List<String>> tableRowData = Lists.newArrayList();
        List<Boolean> shadeVariantData = Lists.newArrayList();

        for(PurpleChrArmCopyNumber chrArmCopyNumber : sort(chrArmCopyNumbers))
        {
            List<String> rowData = Lists.newArrayList();

            tableRowData.add(rowData);

            rowData.add(chrArmCopyNumber.chromosome());
            rowData.add(chrArmCopyNumber.arm());
            rowData.add(chrArmCopyNumber.type());
            rowData.add(formatSingleDigitDecimal(chrArmCopyNumber.copyNumber()));
            rowData.add(formatSingleDigitDecimal(chrArmCopyNumber.relativeCopyNumber()));
            rowData.add(chrArmCopyNumber.driverInterpretation().toString());

            boolean isCandidate = chrArmCopyNumber.driverInterpretation() == DriverInterpretation.LOW;
            shadeVariantData.add(isCandidate);
        }

        int tableWidth = MAX_TABLE_WIDTH - 100;
        scaleTableWidths(headers, tableRowData, tableWidth, reportResources);

        int columnCount = headers.size();
        String[] columnTexts = new String[columnCount];
        float[] columnWidths = new float[columnCount];
        float[] columnPercentages = new float[columnCount];

        convertWidths(headers, columnTexts, columnWidths, columnPercentages);

        BaseTable table = createStandardTable(docCtx, title, tableWidth, columnWidths, columnTexts, reportResources);

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

    private static List<PurpleChrArmCopyNumber> sort(final List<PurpleChrArmCopyNumber> arms)
    {
        return arms.stream().sorted((arm1, arm2) ->
        {
            DriverInterpretation driver1 = arm1.driverInterpretation();
            DriverInterpretation driver2 = arm2.driverInterpretation();

            if(driver1 != driver2)
            {
                return driver1 == DriverInterpretation.HIGH ? -1 : 1;
            }

            int chr1 = HumanChromosome.chromosomeRank(arm1.chromosome());
            int chr2 = HumanChromosome.chromosomeRank(arm2.chromosome());

            if(chr1 != chr2)
            {
                return Integer.compare(chr1, chr2);
            }

            return arm1.arm().compareTo(arm2.arm());

        }).collect(Collectors.toList());
    }
}
