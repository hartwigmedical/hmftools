package com.hartwig.hmftools.orange.report.tables;

import static com.hartwig.hmftools.common.peach.PeachUtil.UNKNOWN_ALLELE_STRING;
import static com.hartwig.hmftools.common.peach.PeachUtil.convertToZygosityString;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.datamodel.peach.PeachGenotype;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class PharmacogeneticsTable
{
    @NotNull
    public static Table build(@NotNull String title, float width, @NotNull Set<PeachGenotype> genotypes,
            @NotNull ReportResources reportResources)
    {
        if(genotypes.isEmpty())
        {
            return new Tables(reportResources).createEmpty(title, width);
        }

        Cells cells = new Cells(reportResources);
        Table contentTable = Tables.createContent(width,
                new float[] { 1, 1, 1, 1, 2, 1 },
                new Cell[] { cells.createHeader("Gene"), cells.createHeader("Haplotype"), cells.createHeader("Genotype"),
                        cells.createHeader("Function"), cells.createHeader("Linked drugs"), cells.createHeader("Source") });

        for(PeachGenotype genotype : sort(genotypes))
        {
            String genotypeString = genotype.allele().equals(UNKNOWN_ALLELE_STRING) ? "" : convertToZygosityString(genotype.alleleCount());

            contentTable.addCell(cells.createContent(genotype.gene()));
            contentTable.addCell(cells.createContent(genotype.allele()));
            contentTable.addCell(cells.createContent(genotypeString));
            contentTable.addCell(cells.createContent(genotype.function()));
            contentTable.addCell(cells.createContent(genotype.linkedDrugs()));
            contentTable.addCell(cells.createUrl(sourceName(genotype.urlPrescriptionInfo()), url(genotype.urlPrescriptionInfo())));
        }

        return new Tables(reportResources).createWrapping(contentTable, title);
    }

    @NotNull
    private static List<PeachGenotype> sort(@NotNull Set<PeachGenotype> genotypes)
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

    @NotNull
    private static String sourceName(@NotNull String urlPrescriptionInfo)
    {
        String url = extractUrl(urlPrescriptionInfo);
        if(url.startsWith("https://www.pharmgkb.org"))
        {
            return "PHARMGKB";
        }
        else
        {
            return Strings.EMPTY;
        }
    }

    @NotNull
    private static String url(@NotNull String urlPrescriptionInfo)
    {
        String url = extractUrl(urlPrescriptionInfo);
        if(url.startsWith("https://www.pharmgkb.org"))
        {
            return url;
        }
        else
        {
            return Strings.EMPTY;
        }
    }

    @NotNull
    private static String extractUrl(@NotNull String urlPrescriptionInfo)
    {
        return urlPrescriptionInfo.split(";")[0];
    }
}
