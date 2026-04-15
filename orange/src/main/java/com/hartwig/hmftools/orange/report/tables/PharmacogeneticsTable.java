package com.hartwig.hmftools.orange.report.tables;

import static com.hartwig.hmftools.common.peach.PeachUtil.UNKNOWN_ALLELE_STRING;
import static com.hartwig.hmftools.common.peach.PeachUtil.convertToZygosityString;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.datamodel.peach.PeachGenotype;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;

import be.quodlibet.boxable.BaseTable;

import com.hartwig.hmftools.orange.report.DocumentContext;

import java.io.IOException;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class PharmacogeneticsTable
{
    public static BaseTable build(final DocumentContext docCtx, final String title, float width, final Set<PeachGenotype> genotypes,
            final ReportResources reportResources) throws IOException
    {
        if(genotypes.isEmpty())
        {
            return new Tables(reportResources).createEmpty(docCtx, title, width);
        }

        float[] colWidths = { 1, 1, 1, 1, 2, 1 };
        String[] headerTexts = { "Gene", "Haplotype", "Genotype", "Function", "Linked drugs", "Source" };

        Cells cells = new Cells(reportResources);
        BaseTable table = TableCommon.createStandardTable(docCtx, title, width, colWidths, headerTexts, reportResources);
        float[] pcts = TableCommon.toPercentages(colWidths);

        for(PeachGenotype genotype : sort(genotypes))
        {
            List<String> rowValues = new ArrayList<>();
            String genotypeString = genotype.allele().equals(UNKNOWN_ALLELE_STRING) ? "" : convertToZygosityString(genotype.alleleCount());

            rowValues.add(genotype.gene());
            rowValues.add(genotype.allele());
            rowValues.add(genotypeString);
            rowValues.add(genotype.function());
            rowValues.add(genotype.linkedDrugs());
            rowValues.add(sourceName(genotype.urlPrescriptionInfo()));
            cells.addRow(table, pcts, rowValues);
        }

        return table;
    }

    private static List<PeachGenotype> sort(final Set<PeachGenotype> genotypes)
    {
        return genotypes.stream().sorted((genotype1, genotype2) ->
        {
            if(!genotype1.gene().equals(genotype2.gene()))
            {
                return genotype1.gene().compareTo(genotype2.gene());
            }
            else if(!genotype1.allele().equals(genotype2.allele()))
            {
                return genotype1.allele().compareTo(genotype2.allele());
            }
            else
            {
                return Integer.compare(genotype1.alleleCount(), genotype2.alleleCount());
            }
        }).collect(Collectors.toList());
    }

    private static final String PHARMGKB_URL = "https://www.pharmgkb.org";
    private static final String PHARMGKB_NAME = "PHARMGKB";

    private static String sourceName(final String urlPrescriptionInfo)
    {
        String url = extractUrl(urlPrescriptionInfo);
        if(url.startsWith(PHARMGKB_URL))
        {
            return PHARMGKB_NAME;
        }
        else
        {
            return Strings.EMPTY;
        }
    }

    @NotNull
    private static String extractUrl(final String urlPrescriptionInfo)
    {
        return urlPrescriptionInfo.split(";")[0];
    }
}
