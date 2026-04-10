package com.hartwig.hmftools.orange.report.tables;

import static java.lang.String.format;

import static com.hartwig.hmftools.orange.algo.OrangeConstants.isCandidateLikelihood;
import static com.hartwig.hmftools.orange.report.tables.SomaticVariantTable.copyNumberDisplay;
import static com.hartwig.hmftools.orange.report.tables.SomaticVariantTable.hgvsDisplay;
import static com.hartwig.hmftools.orange.report.tables.SomaticVariantTable.hotspotDisplay;
import static com.hartwig.hmftools.orange.report.tables.SomaticVariantTable.locationDisplay;
import static com.hartwig.hmftools.orange.report.tables.SomaticVariantTable.sort;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_COPIES;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_GENE;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_HGVS;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.VALUE_NO;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.VALUE_YES;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatSingleDigitDecimal;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatTwoDigitDecimal;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_AF;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_BIALLELIC;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_HOTSPOT;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_MACN;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_VARIANT;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_VCN;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_CN;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_DP;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_DRIVER;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_ZYGOSITY;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.VALUE_HET;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.VALUE_HOM;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.addEntry;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.cellArray;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.floatArray;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.purple.PurpleGenotypeStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleTranscriptImpact;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

public final class GermlineVariantTable
{
    public static final String COL_CLINVAR = "Clinvar";
    public static final String COL_GNOMAD = "Gnomad";

    public static Table build(final String title, float width, final List<PurpleVariant> variants, final ReportResources reportResources)
    {
        if(variants.isEmpty())
        {
            return new Tables(reportResources).createEmpty(title, width);
        }

        Cells cells = new Cells(reportResources);

        List<Float> widths = Lists.newArrayList();
        List<Cell> cellEntries = Lists.newArrayList();

        addEntry(cells, widths, cellEntries, 1.5, COL_GENE);
        addEntry(cells, widths, cellEntries, 2, COL_VARIANT);
        addEntry(cells, widths, cellEntries, 2, COL_HGVS);
        addEntry(cells, widths, cellEntries, 1.5, COL_ZYGOSITY);
        addEntry(cells, widths, cellEntries, 0.75, COL_AF);
        addEntry(cells, widths, cellEntries, 1, COL_DP);
        addEntry(cells, widths, cellEntries, 1, COL_COPIES);
        // addEntry(cells, widths, cellEntries, 0.9, COL_MACN);
        addEntry(cells, widths, cellEntries, 1.5, COL_HOTSPOT);
        addEntry(cells, widths, cellEntries, 1, COL_BIALLELIC);
        addEntry(cells, widths, cellEntries, 1.5, COL_GNOMAD);
        addEntry(cells, widths, cellEntries, 2, COL_CLINVAR);
        addEntry(cells, widths, cellEntries, 1.25, COL_DRIVER);

        Table table = Tables.createContent(width, floatArray(widths), cellArray(cellEntries));

        for(PurpleVariant variant : sort(variants))
        {
            List<PurpleTranscriptImpact> transcriptImpacts = Lists.newArrayList(variant.canonicalImpact());
            transcriptImpacts.addAll(variant.otherImpacts());

            for(PurpleTranscriptImpact transcriptImpact : transcriptImpacts)
            {
                List<Cell> rowCells = Lists.newArrayList();

                rowCells.add(cells.createContent(variant.gene()));
                rowCells.add(cells.createContent(locationDisplay(variant)));
                rowCells.add(cells.createContent(hgvsDisplay(transcriptImpact)));
                rowCells.add(cells.createContent(simplifiedDisplay(variant.genotypeStatus())));
                rowCells.add(cells.createContent(formatTwoDigitDecimal(variant.tumorDepth().alleleFrequency())));
                rowCells.add(cells.createContent(String.valueOf(variant.tumorDepth().totalReadCount())));
                rowCells.add(cells.createContent(copyNumberDisplay(variant)));
                // rowCells.add(cells.createContent(formatSingleDigitDecimal(variant.minorAlleleCopyNumber())));
                rowCells.add(cells.createContent(hotspotDisplay(variant)));
                rowCells.add(cells.createContent(variant.biallelic() ? VALUE_YES : VALUE_NO));
                rowCells.add(cells.createContent(format("%4.2e", variant.gnomadFrequency())));
                rowCells.add(cells.createContent(variant.clinvarPathogenicity()));

                DriverInterpretation driverInterpretation = DriverInterpretation.interpret(variant.driverLikelihood());
                rowCells.add(cells.createContent(driverInterpretation.toString()));

                if(isCandidateLikelihood(variant.driverLikelihood()))
                {
                    reportResources.shadeCandidateCells(rowCells);
                }

                rowCells.forEach(x -> table.addCell(x));
            }
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
