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
import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatTwoDigitDecimal;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_AF;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_BIALLELIC;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_HOTSPOT;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_VARIANT;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_DP;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_DRIVER;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_ZYGOSITY;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.VALUE_HET;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.VALUE_HOM;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.addEntry;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.stringArray;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.createStandardTable;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.createEmptyTable;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.toPercentages;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.floatArray;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.purple.PurpleGenotypeStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleTranscriptImpact;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.Cells;

import be.quodlibet.boxable.Cell;
import be.quodlibet.boxable.BaseTable;

import org.apache.pdfbox.pdmodel.PDPage;

import com.hartwig.hmftools.orange.report.DocumentContext;

import java.io.IOException;

public final class GermlineVariantTable
{
    public static final String COL_CLINVAR = "Clinvar";
    public static final String COL_GNOMAD = "Gnomad";

    public static BaseTable build(final DocumentContext docCtx, final String title, float width, final List<PurpleVariant> variants,
            final ReportResources reportResources) throws IOException
    {
        if(variants.isEmpty())
        {
            return createEmptyTable(docCtx, title, width, reportResources);
        }

        Cells cells = new Cells(reportResources);

        List<Float> widths = Lists.newArrayList();
        List<String> headers = Lists.newArrayList();

        addEntry(widths, headers, 1.5, COL_GENE);
        addEntry(widths, headers, 2, COL_VARIANT);
        addEntry(widths, headers, 2, COL_HGVS);
        addEntry(widths, headers, 1.5, COL_ZYGOSITY);
        addEntry(widths, headers, 0.75, COL_AF);
        addEntry(widths, headers, 1, COL_DP);
        addEntry(widths, headers, 1, COL_COPIES);
        // addEntry(widths, headers, 0.9, COL_MACN);
        addEntry(widths, headers, 1.5, COL_HOTSPOT);
        addEntry(widths, headers, 1, COL_BIALLELIC);
        addEntry(widths, headers, 1.5, COL_GNOMAD);
        addEntry(widths, headers, 2, COL_CLINVAR);
        addEntry(widths, headers, 1.25, COL_DRIVER);

        BaseTable table = createStandardTable(docCtx, title, width, floatArray(widths), stringArray(headers), reportResources);
        float[] pcts = toPercentages(floatArray(widths));

        for(PurpleVariant variant : sort(variants))
        {
            List<PurpleTranscriptImpact> transcriptImpacts = Lists.newArrayList(variant.canonicalImpact());
            transcriptImpacts.addAll(variant.otherImpacts());

            for(PurpleTranscriptImpact transcriptImpact : transcriptImpacts)
            {
                List<String> rowCells = Lists.newArrayList();

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

                List<Cell<PDPage>> createdCells = cells.addRow(table, pcts, rowCells);

                if(isCandidateLikelihood(variant.driverLikelihood()))
                {
                    reportResources.shadeCandidateCells(createdCells);
                }
            }
        }
        return table;
    }

    private static String simplifiedDisplay(PurpleGenotypeStatus genotypeStatus)
    {
        return switch(genotypeStatus)
        {
            case HOM_REF, HOM_ALT -> VALUE_HOM;
            case HET -> VALUE_HET;
            case UNKNOWN -> "UNKNOWN";
        };
    }
}
