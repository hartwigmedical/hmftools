package com.hartwig.hmftools.peach;

import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import org.jetbrains.annotations.NotNull;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class VariantHaplotypeEvent implements HaplotypeEvent
{
    public static final String EVENT_TYPE_STRING = "VAR";
    public static final int ID_FIELD_COUNT = 5;
    @NotNull
    private final String id;
    @NotNull
    public final Chromosome chromosome;
    public final int position;
    @NotNull
    public final String ref;
    @NotNull
    public final String alt;

    public VariantHaplotypeEvent(
            @NotNull final String id,
            @NotNull final Chromosome chromosome,
            final int position,
            @NotNull final String ref,
            @NotNull final String alt
    )
    {
        if (!id.startsWith(EVENT_TYPE_STRING))
        {
            throw new java.lang.IllegalArgumentException(String.format("Invalid ID for VariantHaplotypeEvent: %s", id));
        }

        this.id = id;
        this.chromosome = chromosome;
        this.position = position;
        this.ref = ref;
        this.alt = alt;
    }

    public static VariantHaplotypeEvent fromId(@NotNull String eventId)
    {
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
        return new VariantHaplotypeEvent(eventId, chromosome, position, ref, alt);
    }

    @NotNull
    public String id()
    {
        return this.id;
    }

    public List<Integer> getCoveredPositions(){
        return IntStream.range(position, position + ref.length()).boxed().collect(Collectors.toList());
    }
}
