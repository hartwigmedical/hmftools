package com.hartwig.hmftools.orange.report.tables;

import static java.lang.String.format;

import static com.hartwig.hmftools.orange.algo.OrangeConstants.isCandidateLikelihood;
import static com.hartwig.hmftools.orange.report.ReportResources.formatPercentageField;
import static com.hartwig.hmftools.orange.report.ReportResources.formatSingleDigitDecimal;
import static com.hartwig.hmftools.orange.report.ReportResources.formatTwoDigitDecimal;
import static com.hartwig.hmftools.orange.report.interpretation.Variants.COL_AF;
import static com.hartwig.hmftools.orange.report.interpretation.Variants.COL_BIALLELIC;
import static com.hartwig.hmftools.orange.report.interpretation.Variants.COL_DP;
import static com.hartwig.hmftools.orange.report.interpretation.Variants.COL_HOTSPOT;
import static com.hartwig.hmftools.orange.report.interpretation.Variants.COL_MACN;
import static com.hartwig.hmftools.orange.report.interpretation.Variants.COL_VARIANT;
import static com.hartwig.hmftools.orange.report.interpretation.Variants.COL_VCN;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_CN;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_DRIVER;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_ZYGOSITY;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.VALUE_HET;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.VALUE_HOM;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.addEntry;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.cellArray;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.floatArray;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.intToFloatArray;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.datamodel.purple.PurpleGenotypeStatus;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.datamodel.VariantEntry;
import com.hartwig.hmftools.orange.report.interpretation.Variants;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

public final class GermlineVariantTable
{
    public static final String COL_CLINVAR = "Clinvar";
    public static final String COL_GNOMAD = "Gnomad";

    public static Table build(final String title, float width, final List<VariantEntry> variants, final ReportResources reportResources)
    {
        if(variants.isEmpty())
        {
            return new Tables(reportResources).createEmpty(title, width);
        }

        Cells cells = new Cells(reportResources);

        List<Float> widths = Lists.newArrayList();
        List<Cell> cellEntries = Lists.newArrayList();

        addEntry(cells, widths, cellEntries, 3, COL_VARIANT);
        addEntry(cells, widths, cellEntries, 1.5, COL_ZYGOSITY);
        addEntry(cells, widths, cellEntries, 1, COL_AF);
        addEntry(cells, widths, cellEntries, 1, COL_DP);
        addEntry(cells, widths, cellEntries, 1, COL_VCN);
        addEntry(cells, widths, cellEntries, 1, COL_CN);
        addEntry(cells, widths, cellEntries, 1, COL_MACN);
        addEntry(cells, widths, cellEntries, 1.5, COL_HOTSPOT);
        addEntry(cells, widths, cellEntries, 1.5, COL_BIALLELIC);
        addEntry(cells, widths, cellEntries, 1.5, COL_GNOMAD);
        addEntry(cells, widths, cellEntries, 2, COL_CLINVAR);
        addEntry(cells, widths, cellEntries, 1, COL_DRIVER);

        Table table = Tables.createContent(width, floatArray(widths), cellArray(cellEntries));

        for(VariantEntry variant : Variants.sort(variants))
        {
            List<Cell> rowCells = Lists.newArrayList();

            rowCells.add(cells.createContent(Variants.variantField(variant)));
            rowCells.add(cells.createContent(simplifiedDisplay(variant.genotypeStatus())));
            rowCells.add(cells.createContent(formatTwoDigitDecimal(variant.vaf())));
            rowCells.add(cells.createContent(String.valueOf(variant.depth())));
            rowCells.add(cells.createContent(formatSingleDigitDecimal(variant.variantCopyNumber())));
            rowCells.add(cells.createContent(formatSingleDigitDecimal(variant.totalCopyNumber())));
            rowCells.add(cells.createContent(formatSingleDigitDecimal(variant.minorAlleleCopyNumber())));
            rowCells.add(cells.createContent(Variants.hotspotField(variant)));
            rowCells.add(cells.createContent(variant.biallelic() ? "Yes" : "No"));
            rowCells.add(cells.createContent(format("%4.2e", variant.gnomadFrequency())));
            rowCells.add(cells.createContent(variant.clinvarPathogenicity()));

            rowCells.add(cells.createContent(formatPercentageField(variant.driverLikelihood())));

            if(isCandidateLikelihood(variant.driverLikelihood()))
            {
                reportResources.shadeCandidateCells(rowCells);
            }

            rowCells.forEach(x -> table.addCell(x));
        }

        return new Tables(reportResources).createWrapping(table, title);
    }

    private static String simplifiedDisplay(PurpleGenotypeStatus genotypeStatus)
    {
        switch(genotypeStatus)
        {
            case HOM_REF:
            case HOM_ALT:
                return VALUE_HOM;
            case HET:
                return VALUE_HET;
            case UNKNOWN:
                return "UNKNOWN";
        }
        throw new IllegalStateException();
    }
}
