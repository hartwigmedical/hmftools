package com.hartwig.hmftools.peach.event;

public class HaplotypeEventFactory
{
    public static HaplotypeEvent fromId(final String eventId)
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
