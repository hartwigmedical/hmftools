package com.hartwig.hmftools.orange.report.tables;

import static java.lang.String.format;

import static com.hartwig.hmftools.orange.algo.OrangeConstants.isCandidateLikelihood;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_POSITION;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatSingleDigitDecimal;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_DRIVER;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_GENE;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_JCN;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_TYPE;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_ZYGOSITY;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.VALUE_HET;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.VALUE_HOM;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.addEntry;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.stringArray;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.createStandardTable;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.createEmptyTable;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.toPercentages;
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

import be.quodlibet.boxable.Cell;
import be.quodlibet.boxable.BaseTable;

import org.apache.pdfbox.pdmodel.PDPage;

import com.hartwig.hmftools.orange.report.DocumentContext;

import java.io.IOException;

import org.jetbrains.annotations.Nullable;

public final class DisruptionTable
{
    private static final String COL_CONTEXT = "Context";
    private static final String COL_UNDISRUPTED_CN = "Undisrupted CN";

    public static BaseTable build(final DocumentContext docCtx,
            final String title, float width, final List<LinxBreakend> breakends, final ReportResources reportResources) throws IOException
    {
        if(breakends.isEmpty())
        {
            return createEmptyTable(docCtx, title, width, reportResources);
        }

        Cells cells = new Cells(reportResources);

        List<Integer> widths = Lists.newArrayList();
        List<String> headers = Lists.newArrayList();

        addEntry(widths, headers, 1, COL_GENE);
        addEntry(widths, headers, 3, COL_POSITION);
        addEntry(widths, headers, 1, COL_ZYGOSITY);
        addEntry(widths, headers, 2, COL_CONTEXT);
        addEntry(widths, headers, 1, COL_TYPE);
        addEntry(widths, headers, 1, COL_JCN);
        addEntry(widths, headers, 2, COL_UNDISRUPTED_CN);
        addEntry(widths, headers, 1, COL_DRIVER);

        BaseTable table = createStandardTable(docCtx, title, width, intToFloatArray(widths), stringArray(headers), reportResources);
        float[] pcts = toPercentages(intToFloatArray(widths));

        List<LinxBreakend> sortedBreakends = sort(breakends);

        for(int i = 0; i < sortedBreakends.size(); ++i)
        {
            LinxBreakend breakend = sortedBreakends.get(i);

            LinxBreakend otherBreakend = null;

            if(i + 1 < sortedBreakends.size() && sortedBreakends.get(i + 1).svId() == breakend.svId())
            {
                otherBreakend = sortedBreakends.get(i + 1);
            }

            List<String> rowCells = Lists.newArrayList();

            rowCells.add(cells.createContent(breakend.gene()));
            rowCells.add(cells.createContent(locationDisplay(breakend, otherBreakend)));
            rowCells.add(cells.createContent(breakend.undisruptedCopyNumber() < 0.5 ? VALUE_HOM : VALUE_HET));
            rowCells.add(cells.createContent(contextDisplay(breakend, otherBreakend)));
            rowCells.add(cells.createContent(breakend.type().toString()));
            rowCells.add(cells.createContent(formatSingleDigitDecimal(breakend.junctionCopyNumber())));
            rowCells.add(cells.createContent(formatSingleDigitDecimal(breakend.undisruptedCopyNumber())));

            DriverInterpretation interpretation = DriverInterpretation.interpret(breakend.driverLikelihood());
            rowCells.add(cells.createContent(interpretation.toString()));

            List<Cell<PDPage>> createdCells = cells.addRow(table, pcts, rowCells);

            if(isCandidateLikelihood(breakend.driverLikelihood()))
            {
                reportResources.shadeCandidateCells(createdCells);
            }

            if(otherBreakend != null)
            {
                ++i;
            }
        }

        return table;
    }

    private static String locationDisplay(final LinxBreakend breakend)
    {
        return format("%s:%d", breakend.chromosome(), breakend.position());
    }

    private static String locationDisplay(final LinxBreakend lower, @Nullable final LinxBreakend upper)
    {
        if(upper == null)
        {
            return locationDisplay(lower);
        }

        return format("%s - %s", locationDisplay(lower), locationDisplay(upper));
    }

    private static String contextDisplay(final LinxBreakend lower, @Nullable final LinxBreakend upper)
    {
        if(upper == null)
        {
            return contextStr(lower, true);
        }

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
        return switch(orientation)
        {
            case UPSTREAM -> "Upstream";
            case DOWNSTREAM -> "Downstream";
        };

    }

    private static List<LinxBreakend> sort(final List<LinxBreakend> breakends)
    {
        return breakends.stream().sorted((breakend1, breakend2) ->
        {
            if(breakend1.svId() == breakend2.svId())
            {
                return Integer.compare(breakend1.id(), breakend2.id());
            }

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
