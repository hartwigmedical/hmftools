package com.hartwig.hmftools.pavereverse.serve;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.pavereverse.ReversePave;

import org.jetbrains.annotations.NotNull;

class StatsForGene
{
    final ReversePave baseSequenceVariantsCalculator;
    @NotNull
    public final String mGene;

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

    StatsForGene(final ReversePave baseSequenceVariantsCalculator, @NotNull final String gene)
    {
        this.baseSequenceVariantsCalculator = baseSequenceVariantsCalculator;
        mGene = gene;
    }

    public void recordResultsForAnnotation(VariantStatus comparison)
    {
        if(comparison.mParseException != null)
        {
            System.out.println(comparison.mParseException.getMessage());
        }
        if(comparison.hasProcessingError())
        {
            Throwable t = comparison.mProcessingError;
            System.out.println(t.getMessage());
        }
        if(comparison.hotspotsSame())
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
            Differences.add(new DifferenceWithTransvar(comparison, baseSequenceVariantsCalculator));
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
                "Gene='" + mGene + '\'' +
                ", NumberNotParsed=" + NumberNotParsed +
                ", NumberWithDifferentHotspots=" + NumberWithDifferentHotspots +
                ", NumberWithSameHotspots=" + NumberWithSameHotspots +
                '}';
    }
}
