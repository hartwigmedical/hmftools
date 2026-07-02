package com.hartwig.hmftools.orange.report.tables;

import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_DRIVER;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_TYPE;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.addEntry;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.stringArray;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.createStandardTable;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.createEmptyTable;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.toPercentages;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.floatArray;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatPercentage;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatSingleDigitDecimal;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterEntry;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.Cells;

import be.quodlibet.boxable.Cell;
import be.quodlibet.boxable.BaseTable;

import org.apache.pdfbox.pdmodel.PDPage;

import com.hartwig.hmftools.orange.report.DocumentContext;

import java.io.IOException;

import org.apache.logging.log4j.util.Strings;

public final class ViralPresenceTable
{
    private static final String COL_INTEGRATIONS = "Integrations";
    private static final String COL_QC = "QC Status";
    private static final String COL_PERC_COV = "% Covered";
    private static final String COL_MEAN_COV = "Mean Cov";
    private static final String COL_CLONAL_COV = "Exp Clonal Cov";
    private static final String COL_VIRUS = "Virus";

    public static BaseTable build(final DocumentContext docCtx,
            final String title, float width, final List<VirusInterpreterEntry> viruses, final ReportResources reportResources)
            throws IOException
    {
        if(viruses.isEmpty())
        {
            return createEmptyTable(docCtx, title, width, reportResources);
        }

        Cells cells = new Cells(reportResources);

        List<Float> widths = Lists.newArrayList();
        List<String> headers = Lists.newArrayList();

        addEntry(widths, headers, 2.5, COL_VIRUS);
        addEntry(widths, headers, 2.5, COL_QC);
        addEntry(widths, headers, 1, COL_TYPE);
        addEntry(widths, headers, 2, COL_INTEGRATIONS);
        addEntry(widths, headers, 1.5, COL_PERC_COV);
        addEntry(widths, headers, 1.5, COL_MEAN_COV);
        addEntry(widths, headers, 2, COL_CLONAL_COV);
        addEntry(widths, headers, 1, COL_DRIVER);

        BaseTable table = createStandardTable(docCtx, title, width, floatArray(widths), stringArray(headers), reportResources);
        float[] pcts = toPercentages(floatArray(widths));

        for(VirusInterpreterEntry virus : viruses)
        {
            List<String> rowCells = Lists.newArrayList();

            rowCells.add(cells.createContent(virus.name()));
            rowCells.add(cells.createContent(virus.qcStatus().toString()));
            rowCells.add(cells.createContent(virus.interpretation() != null ? virus.interpretation().name() : Strings.EMPTY));
            rowCells.add(cells.createContent(String.valueOf(virus.integrations())));
            rowCells.add(cells.createContent(formatPercentage(virus.percentageCovered(), false)));
            rowCells.add(cells.createContent(formatSingleDigitDecimal(virus.meanCoverage())));
            rowCells.add(cells.createContent(expectedClonalCoverageField(virus)));
            rowCells.add(cells.createContent(virus.driverInterpretation().toString()));

            List<Cell<PDPage>> createdCells = cells.addRow(table, pcts, rowCells);

            if(virus.driverInterpretation() == DriverInterpretation.LOW)
            {
                reportResources.shadeCandidateCells(createdCells);
            }
        }

        return table;
    }

    private static String expectedClonalCoverageField(final VirusInterpreterEntry virus)
    {
        Double expectedClonalCoverage = virus.expectedClonalCoverage();
        return expectedClonalCoverage != null ? formatSingleDigitDecimal(expectedClonalCoverage) : Strings.EMPTY;
    }
}