package com.hartwig.hmftools.peach;

import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.jetbrains.annotations.NotNull;

import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class VariantHaplotypeEvent implements HaplotypeEvent
{
    public static final String EVENT_TYPE_STRING = "VAR";
    public static final int ID_FIELD_COUNT = 5;
    @NotNull
    public final Chromosome chromosome;
    public final int position;
    @NotNull
    public final String ref;
    @NotNull
    public final String alt;

    public VariantHaplotypeEvent(
            @NotNull final Chromosome chromosome,
            final int position,
            @NotNull final String ref,
            @NotNull final String alt
    )
    {
        this.chromosome = chromosome;
        this.position = position;
        this.ref = ref;
        this.alt = alt;
    }

    public static VariantHaplotypeEvent fromId(@NotNull String eventId)
    {
        if (!eventId.startsWith(EVENT_TYPE_STRING))
        {
            String error_msg = String.format("Invalid ID for VariantHaplotypeEvent: %s", eventId);
            throw new java.lang.IllegalArgumentException(error_msg);
        }

        String[] splitEventId = eventId.split(HaplotypeEvent.EVENT_ID_DELIMITER);
        if (splitEventId.length != VariantHaplotypeEvent.ID_FIELD_COUNT)
        {
            String error_msg = String.format(
                    "ID '%s' of VariantHaplotypeEvent has incorrect field count: %s instead of %s",
                    eventId, splitEventId.length, VariantHaplotypeEvent.ID_FIELD_COUNT
            );
            throw new java.lang.IllegalArgumentException(error_msg);
        }
        Chromosome chromosome = HumanChromosome.fromString(splitEventId[1]);
        int position = Integer.parseInt(splitEventId[2]);
        String ref = splitEventId[3];
        String alt = splitEventId[4];
        VariantHaplotypeEvent event = new VariantHaplotypeEvent(chromosome, position, ref, alt);
        if (!event.id().equals(eventId))
        {
            String error_msg = String.format(
                    "VariantHaplotypeEvent derived from ID ID '%s' has different ID '%s'",
                    eventId, event.id()
            );
            throw new java.lang.IllegalArgumentException(error_msg);
        }
        return event;
    }

    public static VariantHaplotypeEvent fromVariantContext(VariantContext variantContext)
    {
        Chromosome chromosome = HumanChromosome.fromString(variantContext.getContig());
        int position = variantContext.getStart();
        String ref = variantContext.getReference().getBaseString();
        List<Allele> alts = variantContext.getAlternateAlleles();
        if(alts.size() > 1)
        {
            String error_msg = String.format(
                    "Cannot handle variant with multiple alts: %s:%s%s>?",
                    chromosome, position, ref
            );
            throw new IllegalArgumentException(error_msg);
        }
        String alt = variantContext.getAlternateAlleles().get(0).toString();
        return new VariantHaplotypeEvent(chromosome, position, ref, alt);
    }

    @NotNull
    public String id()
    {
        return new StringJoiner(HaplotypeEvent.EVENT_ID_DELIMITER)
                .add(EVENT_TYPE_STRING)
                .add(String.format("chr%1$s", chromosome))
                .add(Integer.toString(position))
                .add(ref)
                .add(alt)
                .toString();
    }

    public List<Integer> getCoveredPositions()
    {
        return IntStream.range(position, position + ref.length()).boxed().collect(Collectors.toList());
    }
}
