package com.hartwig.hmftools.orange.report.tables;

import static java.lang.String.format;

import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatSingleDigitDecimal;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_DRIVER;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_FUSION;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_JCN;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_RNA;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.addEntry;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.cellArray;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.floatArray;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatSupportField;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.linx.FusionPhasedType;
import com.hartwig.hmftools.datamodel.linx.LinxFusion;
import com.hartwig.hmftools.datamodel.common.AllelicDepth;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import com.google.common.collect.Lists;

public final class DnaFusionTable
{
    public static Table build(
            final String title, float width, final List<LinxFusion> fusions, final ReportResources reportResources)
    {
        if(fusions.isEmpty())
        {
            return new Tables(reportResources).createEmpty(title, width);
        }

        Cells cells = new Cells(reportResources);

        List<Float> widths = Lists.newArrayList();
        List<Cell> cellEntries = Lists.newArrayList();

        boolean hasRna = fusions.stream().anyMatch(x -> x.rnaSupport() != null);

        addEntry(cells, widths, cellEntries, 2, COL_FUSION);
        addEntry(cells, widths, cellEntries, 5, "Junctions");
        addEntry(cells, widths, cellEntries, 1, COL_JCN);
        addEntry(cells, widths, cellEntries, 1, "Phasing");
        addEntry(cells, widths, cellEntries, 2, "Type");

        if(hasRna)
            addEntry(cells, widths, cellEntries, 2, COL_RNA);

        addEntry(cells, widths, cellEntries, 1, COL_DRIVER);

        Table table = Tables.createContent(width, floatArray(widths), cellArray(cellEntries));

        for(LinxFusion fusion : sortLinxFusions(fusions))
        {
            List<Cell> rowCells = Lists.newArrayList();

            rowCells.add(cells.createContent(fusionDisplay(fusion)));
            rowCells.add(cells.createContent(transcriptJunctions(fusion)));
            rowCells.add(cells.createContent(formatSingleDigitDecimal(fusion.junctionCopyNumber())));
            rowCells.add(cells.createContent(display(fusion.phased())));
            rowCells.add(cells.createContent(fusion.reportedType().toString()));

            if(hasRna)
            {
                rowCells.add(cells.createContent(rnaSupportField(fusion.rnaSupport())));
            }

            rowCells.add(cells.createContent(fusion.driverInterpretation().toString()));

            if(fusion.driverInterpretation() == DriverInterpretation.LOW)
            {
                reportResources.shadeCandidateCells(rowCells);
            }

            rowCells.forEach(x -> table.addCell(x));
        }

        return new Tables(reportResources).createWrapping(table, title);
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
        switch(fusionPhasedType)
        {
            case INFRAME:
                return "Inframe";
            case SKIPPED_EXONS:
                return "Skipped exons";
            case OUT_OF_FRAME:
                return "Out of frame";
        }
        throw new IllegalStateException();
    }

    private static String rnaSupportField(final AllelicDepth rnaSupport)
    {
        if(rnaSupport == null)
            return ReportResources.NOT_AVAILABLE;
        else
            return formatSupportField(rnaSupport.alleleReadCount(), rnaSupport.totalReadCount());
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

    private static <T> List<T> max5(final List<T> elements)
    {
        return elements.subList(0, Math.min(5, elements.size()));
    }
}
