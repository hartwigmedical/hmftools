package com.hartwig.hmftools.peach.data_loader;

import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.peach.event.VariantHaplotypeEvent;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

import java.io.IOException;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import org.jetbrains.annotations.NotNull;

public class HaplotypeEventLoader
{
    @NotNull
    public static Map<String, Integer> loadRelevantVariantHaplotypeEvents(@NotNull String vcf, @NotNull String sampleName,
            @NotNull Map<Chromosome, Set<Integer>> relevantVariantPositions)
    {
        try(AbstractFeatureReader<VariantContext, LineIterator> reader = getFeatureReader(vcf, new VCFCodec(), false))
        {
            Map<String, Integer> eventIdToCount = new HashMap<>();
            for(VariantContext variantContext : reader.iterator())
            {
                handleVariantContext(variantContext, sampleName, relevantVariantPositions, eventIdToCount);
            }
            return eventIdToCount;
        }
        catch(IOException e)
        {
            throw new RuntimeException(String.format("failed to read VCF file: %s", vcf), e);
        }
    }

    private static void handleVariantContext(@NotNull VariantContext variantContext, @NotNull String sampleName,
            @NotNull Map<Chromosome, Set<Integer>> relevantVariantPositions, @NotNull Map<String, Integer> eventIdToCount)
    {
        Set<Integer> relevantPositionsInChromosome =
                relevantVariantPositions.getOrDefault(HumanChromosome.fromString(variantContext.getContig()), Collections.emptySet());
        VariantHaplotypeEvent event = VariantHaplotypeEvent.fromVariantContext(variantContext);
        Integer count = getEventCount(variantContext.getGenotype(sampleName).getType(), event.id());

        boolean isRelevantEvent = !variantContext.isFiltered() && event.getCoveredPositions()
                .stream()
                .anyMatch(relevantPositionsInChromosome::contains);

        if(isRelevantEvent)
        {
            if(eventIdToCount.containsKey(event.id()))
            {
                throw new IllegalStateException(String.format("encountered event with ID '%s' more than once in input VCF", event.id()));
            }
            eventIdToCount.put(event.id(), count);
        }
    }

    private static Integer getEventCount(@NotNull GenotypeType genotypeType, @NotNull String eventId)
    {
        switch(genotypeType)
        {
            case NO_CALL:
                return null;
            case HOM_REF:
                return 0;
            case HET:
                return 1;
            case HOM_VAR:
                return 2;
            default:
                String errorMsg =
                        String.format("cannot get occurrence count for event with ID '%s' and genotypeType '%s'", eventId, genotypeType);
                throw new IllegalStateException(errorMsg);
        }
    }
}
