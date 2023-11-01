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
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

public class HaplotypeEventLoader
{
    public static Map<String, Integer> loadRelevantVariantHaplotypeEvents(
            String vcf, String sampleName, Map<Chromosome, Set<Integer>> relevantVariantPositions
    )
    {
        try(
                AbstractFeatureReader<VariantContext, LineIterator> reader = getFeatureReader(
                        vcf, new VCFCodec(), false)
        )
        {
            Map<String, Integer> eventIdToCount = new HashMap<>();
            for(VariantContext variantContext : reader.iterator())
            {
                if(variantContext.isFiltered())
                    continue;
                Chromosome chromosome = HumanChromosome.fromString(variantContext.getContig());
                if(!relevantVariantPositions.containsKey(chromosome))
                    continue;
                VariantHaplotypeEvent event = VariantHaplotypeEvent.fromVariantContext(variantContext);
                Set<Integer> relevantPositionsInChromosome = relevantVariantPositions.get(chromosome);
                if(event.getCoveredPositions().stream().noneMatch(relevantPositionsInChromosome::contains))
                    continue;
                Integer count = getEventCount(variantContext.getGenotype(sampleName).getType(), event.id());
                if(count == 0)
                    continue;

                if(eventIdToCount.containsKey(event.id()))
                {
                    throw new IllegalStateException(
                            String.format("encountered event with ID '%s' more than once in VCF '%s'", event.id(), vcf)
                    );
                }

                eventIdToCount.put(event.id(), count);
            }
            return eventIdToCount;
        }
        catch(IOException e)
        {
            throw new RuntimeException(String.format("failed to read VCF file: %s", vcf), e);
        }
    }

    private static Integer getEventCount(GenotypeType genotypeType, String eventId)
    {
        Integer count = null;
        switch(genotypeType)
        {
            case NO_CALL:
            case HOM_REF:
                count = 0;
                break;
            case HET:
                count = 1;
                break;
            case HOM_VAR:
                count = 2;
                break;
            default:
                String error_msg = String.format(
                        "cannot get occurrence count for event with ID '%s' and genotypeType '%s'",
                        eventId,
                        genotypeType
                );
                throw new IllegalStateException(error_msg);
        }
        return count;
    }
}
