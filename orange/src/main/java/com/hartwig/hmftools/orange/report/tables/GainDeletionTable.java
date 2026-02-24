package com.hartwig.hmftools.orange.report.tables;

import static com.hartwig.hmftools.orange.report.ReportResources.formatFoldChangeField;
import static com.hartwig.hmftools.orange.report.ReportResources.formatPercentileField;
import static com.hartwig.hmftools.orange.report.ReportResources.formatSingleDigitDecimal;
import static com.hartwig.hmftools.orange.report.ReportResources.formatTpmField;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.addEntry;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.cellArray;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.floatArray;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.datamodel.purple.PurpleGainDeletion;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.interpretation.Chromosomes;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.apache.logging.log4j.util.Strings;

public final class GainDeletionTable
{
    public static Table build(
            final String title, float width, final List<PurpleGainDeletion> gainsDels, final ReportResources reportResources,
            boolean hasRna)
    {
        if(gainsDels.isEmpty())
        {
            return new Tables(reportResources).createEmpty(title, width);
        }

        Cells cells = new Cells(reportResources);

        List<Integer> widths = Lists.newArrayList();
        List<Cell> cellEntries = Lists.newArrayList();

        addEntry(cells, widths, cellEntries, 1, "Location");
        addEntry(cells, widths, cellEntries, 1, "Gene");
        addEntry(cells, widths, cellEntries, 1, "Type");
        addEntry(cells, widths, cellEntries, 1, "Range"); // if known
        addEntry(cells, widths, cellEntries, 1, "Min CN");
        addEntry(cells, widths, cellEntries, 1, "Max CN");
        addEntry(cells, widths, cellEntries, 1, "Rel CN");
        addEntry(cells, widths, cellEntries, 1, "Driver");

        if(hasRna)
        {
            addEntry(cells, widths, cellEntries, 1, "TPM");
            addEntry(cells, widths, cellEntries, 1, "Percentile");
            addEntry(cells, widths, cellEntries, 1, "Fold Change");
        }

        Table table = Tables.createContent(width, floatArray(widths), cellArray(cellEntries));

        /*
        Table table = Tables.createContent(width,
                new float[] { 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                new Cell[] { cells.createHeader("Location"), cells.createHeader("Gene"),
                        cells.createHeader("Type"), cells.createHeader("CN"), cells.createHeader("TPM"), cells.createHeader("Perc (Type)"),
                        cells.createHeader("FC (Type)"), cells.createHeader("Perc (DB)"), cells.createHeader("FC (DB)") });
        */

        for(PurpleGainDeletion gainDel : sort(gainsDels))
        {
            table.addCell(cells.createContent(gainDel.chromosome() + gainDel.chromosomeBand()));
            table.addCell(cells.createContent(displayGene(gainDel)));
            table.addCell(cells.createContent(gainDel.driver().type().toString()));
            table.addCell(cells.createContent(gainDel.exonRange()));
            table.addCell(cells.createContent(formatSingleDigitDecimal(gainDel.minCopyNumber())));
            table.addCell(cells.createContent(formatSingleDigitDecimal(gainDel.maxCopyNumber())));
            table.addCell(cells.createContent(formatSingleDigitDecimal(gainDel.relativeCopyNumber())));

            table.addCell(cells.createContent(gainDel.driver().driverInterpretation().toString()));

            if(hasRna)
            {
                table.addCell(cells.createContent(formatTpmField(gainDel.tpm())));
                table.addCell(cells.createContent(formatPercentileField(gainDel.tpmPercentile())));
                table.addCell(cells.createContent(formatFoldChangeField(gainDel.tpmFoldChange())));
            }
        }

        return new Tables(reportResources).createWrapping(table, title);
    }

    private static List<PurpleGainDeletion> sort(final List<PurpleGainDeletion> gainsAndDels)
    {
        return gainsAndDels.stream().sorted((gainDel1, gainDel2) ->
        {
            String location1 = Chromosomes.zeroPrefixed(gainDel1.chromosome() + gainDel1.chromosomeBand());
            String location2 = Chromosomes.zeroPrefixed(gainDel2.chromosome() + gainDel2.chromosomeBand());

            if(location1.equals(location2))
            {
                return gainDel1.gene().compareTo(gainDel2.gene());
            }
            else
            {
                return location1.compareTo(location2);
            }
        }).collect(Collectors.toList());
    }

    private static String displayGene(final PurpleGainDeletion gainDel)
    {
        String addon = Strings.EMPTY;
        if(!gainDel.isCanonical())
        {
            addon = " (alt)";
        }
        return gainDel.gene() + addon;
    }
}
