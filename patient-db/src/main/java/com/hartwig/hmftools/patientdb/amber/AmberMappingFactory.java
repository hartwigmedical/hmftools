package com.hartwig.hmftools.patientdb.amber;

import static com.hartwig.hmftools.patientdb.amber.AmberSample.DO_NOT_MATCH;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class AmberMappingFactory
{
    @NotNull
    public static AmberMapping create(final AmberSample firstSample, final AmberSample secondSample)
    {
        final byte[] entries = firstSample.entries();
        byte[] otherEntries = secondSample.entries();
        if(firstSample.sampleId().equals(secondSample.sampleId()))
        {
            throw new IllegalArgumentException("Matching same sample");
        }

        if(entries.length != otherEntries.length)
        {
            throw new IllegalArgumentException("Unable to match different sized identities");
        }

        int matches = 0;
        int sites = 0;
        for(int i = 0; i < entries.length; i++)
        {
            byte myByte = entries[i];
            byte otherByte = otherEntries[i];

            if(myByte != DO_NOT_MATCH && otherByte != DO_NOT_MATCH)
            {
                sites++;
                matches += (myByte == otherByte ? 1 : 0);
            }
        }

        final List<String> sampleNames = Lists.newArrayList(firstSample.sampleId(), secondSample.sampleId());
        Collections.sort(sampleNames);

        return ImmutableAmberMapping.builder()
                .firstSample(sampleNames.get(0))
                .secondSample(sampleNames.get(1))
                .matches(matches)
                .sites(sites)
                .build();
    }

}
