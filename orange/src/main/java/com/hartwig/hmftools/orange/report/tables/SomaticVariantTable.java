package com.hartwig.hmftools.orange.report.tables;

import static java.lang.String.format;

import static com.hartwig.hmftools.orange.algo.OrangeConstants.isCandidateLikelihood;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_COPIES;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_GENE;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_HGVS;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_POSITION;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.VALUE_NO;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.VALUE_YES;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.floatArray;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatPercentageField;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatTwoDigitDecimal;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_AF;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_BIALLELIC;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_CL;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_DRIVER;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_HOTSPOT;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_SL;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_DP;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_RNA;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.addEntry;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.stringArray;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.createStandardTable;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.createEmptyTable;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.toPercentages;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.datamodel.common.AllelicDepth;
import com.hartwig.hmftools.datamodel.purple.PurpleTranscriptImpact;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.Cells;

import be.quodlibet.boxable.Cell;
import be.quodlibet.boxable.BaseTable;

import org.apache.pdfbox.pdmodel.PDPage;

import com.hartwig.hmftools.orange.report.DocumentContext;

import java.io.IOException;

public final class SomaticVariantTable
{
    public static BaseTable build(final DocumentContext docCtx,
            final String title, float width, final List<PurpleVariant> variants, final ReportResources reportResources,
            boolean tumorOnly, boolean hasRna) throws IOException
    {
        if(variants.isEmpty())
        {
            return createEmptyTable(docCtx, title, width, reportResources);
        }

        Cells cells = new Cells(reportResources);

        List<Float> widths = Lists.newArrayList();
        List<String> headers = Lists.newArrayList();

        addEntry(widths, headers, 1.5, COL_GENE);
        addEntry(widths, headers, 2, COL_POSITION);
        addEntry(widths, headers, 3, COL_HGVS);
        addEntry(widths, headers, 1, COL_AF);
        addEntry(widths, headers, 1, COL_DP);
        addEntry(widths, headers, 1, COL_COPIES);
        // addEntry(widths, headers, 1, COL_MACN);
        addEntry(widths, headers, 1.25, COL_HOTSPOT);
        addEntry(widths, headers, 1, COL_BIALLELIC);
        addEntry(widths, headers, 1.25, COL_CL);

        if(tumorOnly)
        {
            addEntry(widths, headers, 1.25, COL_SL);
        }

        if(hasRna)
        {
            addEntry(widths, headers, 1, COL_RNA);
        }

        addEntry(widths, headers, 1, COL_DRIVER);

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
                rowCells.add(cells.createContent(formatTwoDigitDecimal(variant.tumorDepth().alleleFrequency())));
                rowCells.add(cells.createContent(String.valueOf(variant.tumorDepth().totalReadCount())));
                rowCells.add(cells.createContent(copyNumberDisplay(variant)));
                // rowCells.add(cells.createContent(formatSingleDigitDecimal(variant.minorAlleleCopyNumber())));
                rowCells.add(cells.createContent(hotspotDisplay(variant)));
                rowCells.add(cells.createContent(formatPercentageField(variant.biallelicProbability())));

                rowCells.add(cells.createContent(clonalLikelihood(variant)));

                if(tumorOnly)
                {
                    rowCells.add(cells.createContent(variant.somaticLikelihood().toString()));
                }

                if(hasRna)
                {
                    rowCells.add(cells.createContent(rnaDisplay(variant)));
                }

                rowCells.add(cells.createContent(formatPercentageField(variant.driverLikelihood())));

                List<Cell<PDPage>> createdCells = cells.addRow(table, pcts, rowCells);

                if(isCandidateLikelihood(variant.driverLikelihood()))
                {
                    reportResources.shadeCandidateCells(createdCells);
                }

            }
        }

        return table;
    }

    protected static List<PurpleVariant> sort(final List<PurpleVariant> variants)
    {
        return variants.stream().sorted((variant1, variant2) ->
        {
            int driverCompare = Double.compare(variant2.driverLikelihood(), variant1.driverLikelihood());
            if(driverCompare != 0)
            {
                return driverCompare;
            }

            int geneCompare = variant1.gene().compareTo(variant2.gene());
            if(geneCompare != 0)
            {
                return geneCompare;
            }

            return 0;
        }).collect(Collectors.toList());
    }

    protected static String locationDisplay(final PurpleVariant variant)
    {
        boolean phased = variant.localPhaseSets() != null && !variant.localPhaseSets().isEmpty();
        String locationStr = format("%s:%d", variant.chromosome(), variant.position());
        return phased ? format("%s *", locationStr) : locationStr;
    }

    protected static String hgvsDisplay(final PurpleTranscriptImpact transcriptImpact)
    {
        if(transcriptImpact.hgvsProteinImpact().isEmpty())
        {
            return transcriptImpact.hgvsCodingImpact();
        }

        return format("%s [%s]", transcriptImpact.hgvsCodingImpact(), transcriptImpact.hgvsProteinImpact());
    }

    protected static String copyNumberDisplay(final PurpleVariant variant)
    {
        return format("%.1f of %.1f", variant.variantCopyNumber(), variant.adjustedCopyNumber());
    }

    protected static String clonalLikelihood(final PurpleVariant variant)
    {
        return formatPercentageField(1 - variant.subclonalLikelihood());
    }

    protected static String hotspotDisplay(final PurpleVariant variant)
    {
        switch(variant.hotspot())
        {
            case HOTSPOT:
                return VALUE_YES;
            case NEAR_HOTSPOT:
                return "Near";
            default:
                return VALUE_NO;
        }
    }

    protected static String rnaDisplay(final PurpleVariant variant)
    {
        AllelicDepth rnaDepth = variant.rnaDepth();

        if(rnaDepth == null)
        {
            return ReportResources.NOT_AVAILABLE;
        }

        return format("%d/%d", rnaDepth.alleleReadCount(), rnaDepth.totalReadCount());
    }
}

