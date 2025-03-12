package com.hartwig.hmftools.pavereverse.serve;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.pavereverse.ReversePave;

class StatsForGene
{
    final ReversePave Reverse;
    public final String Gene;
    public int NumberProcessedWithoutError = 0;
    public int NumberNotParsed = 0;
    public int NumberWithDifferentHotspots = 0;
    public int NumberWithSameHotspots = 0;
    public int NumberWithSameHotspotsThatUseNonCanonicalTranscript = 0;
    public int NumberWithDifferentHotspotsThatUseNonCanonicalTranscript = 0;
    public int NumberWithNoMatchingTranscript = 0;
    public int NumberForWhichNoVariantCouldBeCalculated = 0;
    public int NumberWithNoUniqueMatchingTranscript = 0;
    public int NumberWithProcessingError = 0;
    public int NumberWithChangesAcrossMultipleExons = 0;

    public final List<VariantStatus> AnnotationsWithDifferentHotspots = new ArrayList<>();
    public final List<DifferenceWithTransvar> Differences = new ArrayList<>();

    StatsForGene(ReversePave baseSequenceVariantsCalculator, String gene)
    {
        Reverse = baseSequenceVariantsCalculator;
        Gene = gene;
    }

    public void recordResultsForAnnotation(VariantStatus comparison)
    {
        if(comparison.ParsingError != null)
        {
            System.out.println(comparison.ParsingError.getMessage());
        }
        if(comparison.hasProcessingError())
        {
            Throwable t = comparison.ProcessingError;
            System.out.println(t.getMessage());
        }
        if(comparison.haveSameChanges())
        {
            NumberWithSameHotspots++;
            if(comparison.usesNonCanonicalTranscript())
            {
                NumberWithSameHotspotsThatUseNonCanonicalTranscript++;
            }
        }
        else
        {
            AnnotationsWithDifferentHotspots.add(comparison);
            NumberWithDifferentHotspots++;
            Differences.add(new DifferenceWithTransvar(comparison, Reverse));
            if(comparison.usesNonCanonicalTranscript())
            {
                NumberWithDifferentHotspotsThatUseNonCanonicalTranscript++;
            }
        }
    }

    public void recordNotParsed()
    {
        NumberNotParsed++;
    }

    @Override
    public String toString()
    {
        return "StatsForGene{" +
                "Gene='" + Gene + '\'' +
                ", NumberNotParsed=" + NumberNotParsed +
                ", NumberWithDifferentHotspots=" + NumberWithDifferentHotspots +
                ", NumberWithSameHotspots=" + NumberWithSameHotspots +
                '}';
    }
}
