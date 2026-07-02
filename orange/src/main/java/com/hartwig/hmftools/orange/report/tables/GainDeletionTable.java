package com.hartwig.hmftools.orange.report.tables;

import static java.lang.String.format;

import static com.hartwig.hmftools.datamodel.purple.PurpleDriverType.GERMLINE_AMP;
import static com.hartwig.hmftools.datamodel.purple.PurpleDriverType.GERMLINE_DELETION;
import static com.hartwig.hmftools.orange.report.ReportResources.NOT_AVAILABLE;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.VALUE_HET;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.VALUE_HOM;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatFoldChangeField;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatPercentileField;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatSingleDigitDecimal;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatTpmField;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_DRIVER;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_GENE;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_LOCATION;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_RANGE;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_REL_CN;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_TPM;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_TYPE;
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
import com.hartwig.hmftools.datamodel.purple.PurpleGainDeletion;
import com.hartwig.hmftools.datamodel.purple.PurpleGermlineStatus;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.Cells;

import be.quodlibet.boxable.Cell;
import be.quodlibet.boxable.BaseTable;

import org.apache.pdfbox.pdmodel.PDPage;

import com.hartwig.hmftools.orange.report.DocumentContext;

import java.io.IOException;

import org.apache.logging.log4j.util.Strings;

public final class GainDeletionTable
{
    private static String COL_MIN_CN = "Min CN";
    private static String COL_MAX_CN = "Max CN";
    private static String COL_ARM_CN = "Arm CN";

    public static BaseTable build(final DocumentContext docCtx,
            final String title, float width, final List<PurpleGainDeletion> gainsDels, final ReportResources reportResources,
            boolean hasRna) throws IOException
    {
        if(gainsDels.isEmpty())
        {
            return createEmptyTable(docCtx, title, width, reportResources);
        }

        Cells cells = new Cells(reportResources);

        List<Integer> widths = Lists.newArrayList();
        List<String> headers = Lists.newArrayList();

        addEntry(widths, headers, 1, COL_LOCATION);
        addEntry(widths, headers, 1, COL_GENE);
        addEntry(widths, headers, 1, COL_TYPE);
        addEntry(widths, headers, 1, COL_RANGE); // if known
        addEntry(widths, headers, 1, COL_MIN_CN);
        addEntry(widths, headers, 1, COL_MAX_CN);
        addEntry(widths, headers, 1, COL_REL_CN);
        addEntry(widths, headers, 1, COL_ARM_CN);

        boolean anyEventsHaveRna = hasRna && gainsDels.stream().anyMatch(x -> x.tpm() != null);

        if(anyEventsHaveRna)
        {
            addEntry(widths, headers, 1, COL_TPM);
            addEntry(widths, headers, 1, "Percentile");
            addEntry(widths, headers, 1, "Fold Change");
        }

        addEntry(widths, headers, 1, COL_DRIVER);

        BaseTable table = createStandardTable(docCtx, title, width, intToFloatArray(widths), stringArray(headers), reportResources);
        float[] pcts = toPercentages(intToFloatArray(widths));

        for(PurpleGainDeletion gainDel : sort(gainsDels))
        {
            List<String> rowCells = Lists.newArrayList();

            rowCells.add(cells.createContent(gainDel.chromosome() + gainDel.chromosomeBand()));
            rowCells.add(cells.createContent(geneDisplay(gainDel)));
            rowCells.add(cells.createContent(typeDisplay(gainDel)));
            rowCells.add(cells.createContent(rangeDisplay(gainDel)));
            rowCells.add(cells.createContent(formatSingleDigitDecimal(gainDel.minCopyNumber())));
            rowCells.add(cells.createContent(formatSingleDigitDecimal(gainDel.maxCopyNumber())));
            rowCells.add(cells.createContent(formatSingleDigitDecimal(gainDel.relativeCopyNumber())));
            rowCells.add(cells.createContent(formatSingleDigitDecimal(gainDel.armCopyNumber())));

            if(anyEventsHaveRna)
            {
                if(gainDel.tpm() != null)
                {
                    rowCells.add(cells.createContent(formatTpmField(gainDel.tpm())));
                    rowCells.add(cells.createContent(formatPercentileField(gainDel.tpmPercentile())));
                    rowCells.add(cells.createContent(formatFoldChangeField(gainDel.tpmFoldChange())));
                }
                else
                {
                    rowCells.add(cells.createContent(NOT_AVAILABLE));
                    rowCells.add(cells.createContent(NOT_AVAILABLE));
                    rowCells.add(cells.createContent(NOT_AVAILABLE));
                }
            }

            rowCells.add(cells.createContent(gainDel.driver().driverInterpretation().toString()));

            List<Cell<PDPage>> createdCells = cells.addRow(table, pcts, rowCells);

            if(gainDel.driver().driverInterpretation() == DriverInterpretation.LOW)
            {
                reportResources.shadeCandidateCells(createdCells);
            }

        }

        return table;
    }

    private static List<PurpleGainDeletion> sort(final List<PurpleGainDeletion> gainsAndDels)
    {
        return gainsAndDels.stream().sorted((gainDel1, gainDel2) ->
        {
            double likelihood1 = gainDel1.driver().driverLikelihood();
            double likelihood2 = gainDel2.driver().driverLikelihood();

            int likelihoodCompare = Double.compare(-likelihood1, -likelihood2);

            if(likelihoodCompare != 0)
            {
                return likelihoodCompare;
            }

            String location1 = zeroPrefixed(gainDel1.chromosome() + gainDel1.chromosomeBand());
            String location2 = zeroPrefixed(gainDel2.chromosome() + gainDel2.chromosomeBand());

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

    private static String rangeDisplay(final PurpleGainDeletion gainDel)
    {
        if(gainDel.exonStart() == null && gainDel.exonEnd() == null)
        {
            return gainDel.geneRange();
        }

        return format("Exon %d - Exon %d", gainDel.exonStart(), gainDel.exonEnd());
    }

    private static String typeDisplay(final PurpleGainDeletion gainDel)
    {
        if(gainDel.driver().type() == GERMLINE_AMP)
        {
            return "AMP";
        }
        else if(gainDel.driver().type() == GERMLINE_DELETION)
        {
            return gainDel.germlineAmpDelFields().germlineStatus() == PurpleGermlineStatus.HOM_DELETION ? VALUE_HOM : VALUE_HET;
        }

        return gainDel.driver().type().toString();
    }

    private static String geneDisplay(final PurpleGainDeletion gainDel)
    {
        String addon = Strings.EMPTY;
        if(!gainDel.isCanonical())
        {
            addon = " (alt)";
        }
        return gainDel.gene() + addon;
    }
}
