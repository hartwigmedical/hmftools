package com.hartwig.hmftools.orange.report.tables;

import static com.hartwig.hmftools.orange.algo.OrangeConstants.isCandidateLikelihood;
import static com.hartwig.hmftools.orange.report.ReportResources.formatSingleDigitDecimal;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_DRIVER;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_GENE;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_JCN;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_LOCATION;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_RANGE;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_TYPE;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_ZYGOSITY;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.VALUE_HET;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.VALUE_HOM;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.addEntry;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.cellArray;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.intToFloatArray;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.datamodel.BreakendEntry;
import com.hartwig.hmftools.orange.report.interpretation.Chromosomes;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.apache.logging.log4j.util.Strings;

public final class DisruptionTable
{
    public static Table build(
            final String title, float width, final List<BreakendEntry> breakends, final ReportResources reportResources)
    {
        if(breakends.isEmpty())
        {
            return new Tables(reportResources).createEmpty(title, width);
        }

        Cells cells = new Cells(reportResources);

        List<Integer> widths = Lists.newArrayList();
        List<Cell> cellEntries = Lists.newArrayList();

        addEntry(cells, widths, cellEntries, 1, COL_LOCATION);
        addEntry(cells, widths, cellEntries, 1, COL_GENE);
        addEntry(cells, widths, cellEntries, 1, COL_ZYGOSITY);
        addEntry(cells, widths, cellEntries, 2, COL_RANGE);
        addEntry(cells, widths, cellEntries, 1, COL_TYPE);
        addEntry(cells, widths, cellEntries, 1, COL_JCN);
        addEntry(cells, widths, cellEntries, 2, "Undisrupted CN");
        addEntry(cells, widths, cellEntries, 1, COL_DRIVER);

        Table table = Tables.createContent(width, intToFloatArray(widths), cellArray(cellEntries));

        for(BreakendEntry breakend : sort(breakends))
        {
            List<Cell> rowCells = Lists.newArrayList();

            rowCells.add(cells.createContent(breakend.location()));
            rowCells.add(cells.createContent(displayGene(breakend)));
            rowCells.add(cells.createContent(breakend.undisruptedCopyNumber() < 0.5 ? VALUE_HOM : VALUE_HET));
            rowCells.add(cells.createContent(breakend.range()));
            rowCells.add(cells.createContent(breakend.type().toString()));
            rowCells.add(cells.createContent(formatSingleDigitDecimal(breakend.junctionCopyNumber())));
            rowCells.add(cells.createContent(formatSingleDigitDecimal(breakend.undisruptedCopyNumber())));

            DriverInterpretation interpretation = DriverInterpretation.interpret(breakend.driverLikelihood());
            rowCells.add(cells.createContent(interpretation.toString()));

            if(isCandidateLikelihood(breakend.driverLikelihood()))
            {
                reportResources.shadeCandidateCells(rowCells);
            }

            rowCells.forEach(x -> table.addCell(x));
        }

        return new Tables(reportResources).createWrapping(table, title);
    }

    private static List<BreakendEntry> sort(final List<BreakendEntry> breakends)
    {
        return breakends.stream().sorted((breakend1, breakend2) ->
        {
            String location1 = Chromosomes.zeroPrefixed(breakend1.location());
            String location2 = Chromosomes.zeroPrefixed(breakend2.location());

            int locationCompare = location1.compareTo(location2);
            if(locationCompare != 0)
            {
                return locationCompare;
            }

            int geneCompare = breakend1.gene().compareTo(breakend2.gene());
            if(geneCompare != 0)
            {
                return geneCompare;
            }

            return breakend1.exonUp() - breakend2.exonUp();
        }).collect(Collectors.toList());
    }

    private static String displayGene(final BreakendEntry breakend)
    {
        String addon = Strings.EMPTY;
        if(!breakend.canonical())
        {
            addon = " (alt)";
        }
        return breakend.gene() + addon;
    }
}
