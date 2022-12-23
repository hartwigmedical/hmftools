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

import static com.hartwig.hmftools.peach.PeachUtils.PCH_LOGGER;
import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

public class HaplotypeEventLoader
{
    public static Map<String, Integer> loadRelevantVariantHaplotypeEvents(
            String vcf, String sampleName, Map<Chromosome, Set<Integer>> relevantVariantPositions
    )
    {
        Map<String, Integer> eventIdToCount = new HashMap<>();
        try(
                AbstractFeatureReader<VariantContext, LineIterator> reader = getFeatureReader(
                        vcf, new VCFCodec(), false)
        )
        {
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
                    PCH_LOGGER.error("encountered event with ID '{}' more than once in VCF '{}'", event.id(), vcf);
                    System.exit(1);
                }

                eventIdToCount.put(event.id(), count);
            }
        }
        catch(IOException e)
        {
            PCH_LOGGER.error("failed to read VCF file({}): {}", vcf, e.toString());
            System.exit(1);
        }
        return eventIdToCount;
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
                PCH_LOGGER.error("cannot get occurrence count for event with ID '{}' with genotypeType '{}'", eventId, genotypeType.toString());
                System.exit(1);
        }
        return count;
    }
}
