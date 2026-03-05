package com.hartwig.hmftools.orange.report.tables;

import static com.hartwig.hmftools.orange.report.ReportResources.formatSingleDigitDecimal;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_CHR;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_CN;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_REL_CN;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_TYPE;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.addEntry;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.cellArray;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.intToFloatArray;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.datamodel.purple.PurpleChrArmCopyNumber;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

public final class ChrArmCopyNumberTable
{
    public static Table build(
            final String title, float width, final List<PurpleChrArmCopyNumber> chrArmCopyNumbers, final ReportResources reportResources)
    {
        if(chrArmCopyNumbers.isEmpty())
        {
            return new Tables(reportResources).createEmpty(title, width);
        }

        Cells cells = new Cells(reportResources);

        List<Integer> widths = Lists.newArrayList();
        List<Cell> cellEntries = Lists.newArrayList();

        addEntry(cells, widths, cellEntries, 1, COL_CHR);
        addEntry(cells, widths, cellEntries, 1, "Arm");
        addEntry(cells, widths, cellEntries, 1, COL_TYPE);
        addEntry(cells, widths, cellEntries, 1, COL_CN);
        addEntry(cells, widths, cellEntries, 1, COL_REL_CN);

        Table table = Tables.createContent(width, intToFloatArray(widths), cellArray(cellEntries));

        for(PurpleChrArmCopyNumber chrArmCopyNumber : chrArmCopyNumbers)
        {
            table.addCell(cells.createContent(chrArmCopyNumber.chromosome()));
            table.addCell(cells.createContent(chrArmCopyNumber.arm()));
            table.addCell(cells.createContent(chrArmCopyNumber.type()));
            table.addCell(cells.createContent(formatSingleDigitDecimal(chrArmCopyNumber.copyNumber())));
            table.addCell(cells.createContent(formatSingleDigitDecimal(chrArmCopyNumber.relativeCopyNumber())));
        }

        return new Tables(reportResources).createWrapping(table, title);
    }
}
