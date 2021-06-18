package com.hartwig.hmftools.neo.cohort;

import static com.hartwig.hmftools.neo.cohort.AlleleCoverage.EXPECTED_ALLELE_COUNT;
import static com.hartwig.hmftools.neo.cohort.AlleleCoverage.gene;

import java.util.List;

public class PredictionData
{
    public final String Allele;
    public final int NeId;
    public final String Peptide;
    public final double Affinity;
    public final double Presentation;

    public static final String DELIM = ",";

    public PredictionData(final String allele, int neId, String peptide, double affinity, double presentation)
    {
        Allele = allele;
        NeId = neId;
        Peptide = peptide;
        Affinity = affinity;
        Presentation = presentation;
    }

    public String toString()
    {
        return String.format("allele(%s) neId(%d) pep(%s) affinity(%.1f) pres(%.4f)",
                Allele, NeId, Peptide, Affinity, Presentation);
    }

    public static PredictionData fromCsv(
            final String data, int alleleIndex, int neIdIndex, int peptideIndex, int affinityIndex, int presentationIndex)
    {
        final String[] items = data.split(DELIM, -1);

        String allele = items[alleleIndex].replaceAll("HLA-", "");

        return new PredictionData(
                allele, Integer.parseInt(items[neIdIndex]), items[peptideIndex],
                Double.parseDouble(items[affinityIndex]), Double.parseDouble(items[presentationIndex]));
    }

    public static void expandHomozygous(final List<PredictionData> predictions)
    {
        if(predictions.size() >= EXPECTED_ALLELE_COUNT)
            return;

        int index = 0;

        while(index < predictions.size())
        {
            PredictionData coverage = predictions.get(index);
            PredictionData nextCoverage = index < predictions.size() - 1 ? predictions.get(index + 1) : null;

            if(nextCoverage != null && gene(coverage.Allele).equals(gene(nextCoverage.Allele)))
            {
                // both alleles present
                index += 2;
                continue;
            }

            if(nextCoverage != null)
            {
                // must be A or B
                predictions.add(index, coverage); // replicate
                index += 2;
                continue;
            }
            else
            {
                predictions.add(coverage);
                break;
            }
        }
    }

}
