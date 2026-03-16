package com.hartwig.hmftools.orange.report.tables;

import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_DRIVER;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.addEntry;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.cellArray;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatPercentage;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatSingleDigitDecimal;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.intToFloatArray;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterEntry;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.apache.logging.log4j.util.Strings;

public final class ViralPresenceTable
{
    public static Table build(
            final String title, float width, final List<VirusInterpreterEntry> viruses, final ReportResources reportResources)
    {
        if(viruses.isEmpty())
        {
            return new Tables(reportResources).createEmpty(title, width);
        }

        Cells cells = new Cells(reportResources);

        List<Integer> widths = Lists.newArrayList();
        List<Cell> cellEntries = Lists.newArrayList();

        addEntry(cells, widths, cellEntries, 3, "Virus");
        addEntry(cells, widths, cellEntries, 3, "QC Status");
        addEntry(cells, widths, cellEntries, 1, "Type");
        addEntry(cells, widths, cellEntries, 1, "Int");
        addEntry(cells, widths, cellEntries, 2, "% Covered");
        addEntry(cells, widths, cellEntries, 2, "Mean Cov");
        addEntry(cells, widths, cellEntries, 2, "Exp Clon Cov");
        addEntry(cells, widths, cellEntries, 2, COL_DRIVER);

        Table table = Tables.createContent(width, intToFloatArray(widths), cellArray(cellEntries));

        for(VirusInterpreterEntry virus : viruses)
        {
            List<Cell> rowCells = Lists.newArrayList();

            rowCells.add(cells.createContent(virus.name()));
            rowCells.add(cells.createContent(virus.qcStatus().toString()));
            rowCells.add(cells.createContent(virus.interpretation() != null ? virus.interpretation().name() : Strings.EMPTY));
            rowCells.add(cells.createContent(String.valueOf(virus.integrations())));
            rowCells.add(cells.createContent(formatPercentage(virus.percentageCovered(), false)));
            rowCells.add(cells.createContent(formatSingleDigitDecimal(virus.meanCoverage())));
            rowCells.add(cells.createContent(expectedClonalCoverageField(virus)));
            rowCells.add(cells.createContent(virus.driverInterpretation().toString()));

            if(virus.driverInterpretation() == DriverInterpretation.LOW)
            {
                reportResources.shadeCandidateCells(rowCells);
            }

            rowCells.forEach(x -> table.addCell(x));
        }

        return new Tables(reportResources).createWrapping(table, title);
    }

    private static String expectedClonalCoverageField(final VirusInterpreterEntry virus)
    {
        Double expectedClonalCoverage = virus.expectedClonalCoverage();
        return expectedClonalCoverage != null ? formatSingleDigitDecimal(expectedClonalCoverage) : Strings.EMPTY;
    }
}