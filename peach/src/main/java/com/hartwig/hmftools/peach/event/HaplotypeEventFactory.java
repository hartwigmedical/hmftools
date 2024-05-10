package com.hartwig.hmftools.peach.event;

import org.jetbrains.annotations.NotNull;

public class HaplotypeEventFactory
{
    @NotNull
    public static HaplotypeEvent fromId(@NotNull String eventId)
    {
        String[] splitEventId = eventId.split(HaplotypeEvent.EVENT_ID_DELIMITER);
        String eventTypeId = splitEventId[0];
        if(eventTypeId.equals(VariantHaplotypeEvent.EVENT_TYPE_STRING))
        {
            return VariantHaplotypeEvent.fromId(eventId);
        }
        else
        {
            String errorMsg = String.format("Cannot construct HaplotypeEvent from id: %s", eventId);
            throw new java.lang.IllegalArgumentException(errorMsg);
        }
    }
}
