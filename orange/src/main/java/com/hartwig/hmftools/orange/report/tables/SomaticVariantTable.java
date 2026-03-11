package com.hartwig.hmftools.orange.report.tables;

import static com.hartwig.hmftools.orange.algo.OrangeConstants.isCandidateLikelihood;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatPercentageField;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatSingleDigitDecimal;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatTwoDigitDecimal;
import static com.hartwig.hmftools.orange.report.interpretation.Variants.COL_AF;
import static com.hartwig.hmftools.orange.report.interpretation.Variants.COL_BIALLELIC;
import static com.hartwig.hmftools.orange.report.interpretation.Variants.COL_CL;
import static com.hartwig.hmftools.orange.report.interpretation.Variants.COL_DRIVER;
import static com.hartwig.hmftools.orange.report.interpretation.Variants.COL_HOTSPOT;
import static com.hartwig.hmftools.orange.report.interpretation.Variants.COL_MACN;
import static com.hartwig.hmftools.orange.report.interpretation.Variants.COL_SL;
import static com.hartwig.hmftools.orange.report.interpretation.Variants.COL_VARIANT;
import static com.hartwig.hmftools.orange.report.interpretation.Variants.COL_VCN;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_CN;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_DP;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_RNA;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.addEntry;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.cellArray;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.intToFloatArray;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.datamodel.VariantEntry;
import com.hartwig.hmftools.orange.report.interpretation.Variants;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

public final class SomaticVariantTable
{
    public static Table build(
            final String title, float width, final List<VariantEntry> variants, final ReportResources reportResources,
            boolean tumorOnly, boolean hasRna)
    {
        if(variants.isEmpty())
        {
            return new Tables(reportResources).createEmpty(title, width);
        }

        Cells cells = new Cells(reportResources);

        List<Integer> widths = Lists.newArrayList();
        List<Cell> cellEntries = Lists.newArrayList();

        addEntry(cells, widths, cellEntries, 3, COL_VARIANT);
        addEntry(cells, widths, cellEntries, 1, COL_AF);
        addEntry(cells, widths, cellEntries, 1, COL_DP);
        addEntry(cells, widths, cellEntries, 1, COL_VCN);
        addEntry(cells, widths, cellEntries, 1, COL_CN);
        addEntry(cells, widths, cellEntries, 1, COL_MACN);
        addEntry(cells, widths, cellEntries, 1, COL_HOTSPOT);
        addEntry(cells, widths, cellEntries, 1, COL_BIALLELIC);
        addEntry(cells, widths, cellEntries, 1, COL_CL);

        if(tumorOnly)
        {
            addEntry(cells, widths, cellEntries, 1, COL_SL);
        }

        if(hasRna)
        {
            addEntry(cells, widths, cellEntries, 1, COL_RNA);
        }

        addEntry(cells, widths, cellEntries, 1, COL_DRIVER);

        Table table = Tables.createContent(width, intToFloatArray(widths), cellArray(cellEntries));

        for(VariantEntry variant : Variants.sort(variants))
        {
            List<Cell> rowCells = Lists.newArrayList();

            rowCells.add(cells.createContent(Variants.variantField(variant)));
            rowCells.add(cells.createContent(formatTwoDigitDecimal(variant.vaf())));
            rowCells.add(cells.createContent(String.valueOf(variant.depth())));
            rowCells.add(cells.createContent(formatSingleDigitDecimal(variant.variantCopyNumber())));
            rowCells.add(cells.createContent(formatSingleDigitDecimal(variant.totalCopyNumber())));
            rowCells.add(cells.createContent(formatSingleDigitDecimal(variant.minorAlleleCopyNumber())));
            rowCells.add(cells.createContent(Variants.hotspotField(variant)));
            rowCells.add(cells.createContent(formatPercentageField(variant.biallelicProbability())));

            rowCells.add(cells.createContent(formatPercentageField(variant.clonalLikelihood())));

            if(tumorOnly)
                rowCells.add(cells.createContent(variant.somaticLikelihood()));

            if(hasRna)
                rowCells.add(cells.createContent(Variants.rnaInfoField(variant)));

            rowCells.add(cells.createContent(formatPercentageField(variant.driverLikelihood())));

            if(isCandidateLikelihood(variant.driverLikelihood()))
            {
                reportResources.shadeCandidateCells(rowCells);
            }

            rowCells.forEach(x -> table.addCell(x));
        }

        return new Tables(reportResources).createWrapping(table, title);
    }
}

