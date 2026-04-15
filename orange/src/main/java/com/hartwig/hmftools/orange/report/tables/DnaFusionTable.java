package com.hartwig.hmftools.orange.report.tables;

import static java.lang.String.format;

import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_JUNCTIONS;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_RNA_FRAGS;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_TYPE;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatSingleDigitDecimal;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_DRIVER;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_FUSION;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_JCN;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.addEntry;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.stringArray;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.createStandardTable;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.createEmptyTable;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.toPercentages;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.floatArray;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.linx.FusionPhasedType;
import com.hartwig.hmftools.datamodel.linx.LinxFusion;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.Cells;

import be.quodlibet.boxable.Cell;
import be.quodlibet.boxable.BaseTable;

import org.apache.pdfbox.pdmodel.PDPage;

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

        Cells cells = new Cells(reportResources);

        List<Float> widths = Lists.newArrayList();
        List<String> headers = Lists.newArrayList();

        boolean hasRna = fusions.stream().anyMatch(x -> x.rnaSupport() != null);

        addEntry(widths, headers, 2, COL_FUSION);
        addEntry(widths, headers, 5, COL_JUNCTIONS);
        addEntry(widths, headers, 1, COL_JCN);
        addEntry(widths, headers, 1, COL_PHASING);
        addEntry(widths, headers, 3, COL_TYPE);

        if(hasRna)
        {
            addEntry(widths, headers, 1, COL_RNA_FRAGS);
        }

        addEntry(widths, headers, 1, COL_DRIVER);

        BaseTable table = createStandardTable(docCtx, title, width, floatArray(widths), stringArray(headers), reportResources);
        float[] pcts = toPercentages(floatArray(widths));

        for(LinxFusion fusion : sortLinxFusions(fusions))
        {
            List<String> rowCells = Lists.newArrayList();

            rowCells.add(cells.createContent(fusionDisplay(fusion)));
            rowCells.add(cells.createContent(transcriptJunctions(fusion)));
            rowCells.add(cells.createContent(formatSingleDigitDecimal(fusion.junctionCopyNumber())));
            rowCells.add(cells.createContent(display(fusion.phased())));
            rowCells.add(cells.createContent(fusion.reportedType().toString()));

            if(hasRna)
            {
                rowCells.add(cells.createContent(String.valueOf(fusion.rnaSupport().alleleReadCount())));
            }

            rowCells.add(cells.createContent(fusion.driverInterpretation().toString()));

            List<Cell<PDPage>> createdCells = cells.addRow(table, pcts, rowCells);

            if(fusion.driverInterpretation() == DriverInterpretation.LOW)
            {
                reportResources.shadeCandidateCells(createdCells);
            }

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

    /*
    private static String rnaSupportField(final AllelicDepth rnaSupport)
    {
        if(rnaSupport == null)
            return ReportResources.NOT_AVAILABLE;
        else
            return formatSupportField(rnaSupport.alleleReadCount(), rnaSupport.totalReadCount());
    }
    */

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
