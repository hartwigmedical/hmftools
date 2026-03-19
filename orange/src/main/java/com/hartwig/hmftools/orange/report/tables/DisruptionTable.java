package com.hartwig.hmftools.orange.report.tables;

import static java.lang.String.format;

import static com.hartwig.hmftools.orange.algo.OrangeConstants.isCandidateLikelihood;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatSingleDigitDecimal;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_DRIVER;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_GENE;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_JCN;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_LOCATION;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_TYPE;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_ZYGOSITY;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.VALUE_HET;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.VALUE_HOM;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.addEntry;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.cellArray;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.intToFloatArray;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.zeroPrefixed;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.linx.LinxBreakend;
import com.hartwig.hmftools.datamodel.linx.LinxGeneOrientation;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.jetbrains.annotations.Nullable;

public final class DisruptionTable
{
    private static final String COL_CONTEXT = "Context";
    private static final String COL_UNDISRUPTED_CN = "Undisrupted CN";

    public static Table build(
            final String title, float width, final List<LinxBreakend> breakends, final ReportResources reportResources)
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
        addEntry(cells, widths, cellEntries, 2, COL_CONTEXT);
        addEntry(cells, widths, cellEntries, 1, COL_TYPE);
        addEntry(cells, widths, cellEntries, 1, COL_JCN);
        addEntry(cells, widths, cellEntries, 2, COL_UNDISRUPTED_CN);
        addEntry(cells, widths, cellEntries, 1, COL_DRIVER);

        Table table = Tables.createContent(width, intToFloatArray(widths), cellArray(cellEntries));

        List<LinxBreakend> sortedBreakends = sort(breakends);

        for(int i = 0; i < sortedBreakends.size(); ++i)
        {
            LinxBreakend breakend = sortedBreakends.get(i);

            LinxBreakend otherBreakend = null;

            if(i + 1 < sortedBreakends.size() && sortedBreakends.get(i + 1).svId() == breakend.svId())
            {
                otherBreakend = sortedBreakends.get(i + 1);
            }

            List<Cell> rowCells = Lists.newArrayList();

            rowCells.add(cells.createContent(locationDisplay(breakend)));
            rowCells.add(cells.createContent(breakend.gene()));
            rowCells.add(cells.createContent(breakend.undisruptedCopyNumber() < 0.5 ? VALUE_HOM : VALUE_HET));
            rowCells.add(cells.createContent(contextDisplay(breakend, otherBreakend)));
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

            if(otherBreakend != null)
                ++i;
        }

        return new Tables(reportResources).createWrapping(table, title);
    }

    private static String locationDisplay(final LinxBreakend breakend)
    {
        return format("%s %s", breakend.chromosome(), breakend.chromosomeBand());
    }

    private static String contextDisplay(final LinxBreakend lower, @Nullable final LinxBreakend upper)
    {
        if(upper == null)
            return contextStr(lower, true);

        return format("%s - %s", contextStr(lower, false), contextStr(upper, false));
    }

    private static String contextStr(final LinxBreakend breakend, boolean includeOrientation)
    {
        String exonRange = null;
        if(breakend.exonUp() > 0)
        {
            if(breakend.exonUp() == breakend.exonDown())
            {
                exonRange = String.format("Exon %d", breakend.exonUp());
            }
            else if(breakend.exonDown() - breakend.exonUp() == 1)
            {
                exonRange = String.format("Intron %d", breakend.exonUp());
            }
        }
        else if(breakend.exonUp() == 0 && (breakend.exonDown() == 1 || breakend.exonDown() == 2))
        {
            exonRange = "Promoter Region";
        }

        return includeOrientation ? format("%s %s", exonRange, orientationStr(breakend.geneOrientation())) : exonRange;
    }

    private static String orientationStr(final LinxGeneOrientation orientation)
    {
        switch(orientation)
        {
            case UPSTREAM: return "Upstream";
            case DOWNSTREAM: return "Downstream";
        }

        return null;
    }

    private static List<LinxBreakend> sort(final List<LinxBreakend> breakends)
    {
        return breakends.stream().sorted((breakend1, breakend2) ->
        {
            if(breakend1.svId() == breakend2.svId())
                return Integer.compare(breakend1.id(), breakend2.id());

            String location1 = zeroPrefixed(breakend1.chromosome()) + breakend1.chromosomeBand();
            String location2 = zeroPrefixed(breakend2.chromosome()) + breakend2.chromosomeBand();

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
}
