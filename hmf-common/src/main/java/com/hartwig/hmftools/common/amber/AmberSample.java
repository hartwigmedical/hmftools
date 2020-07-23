package com.hartwig.hmftools.common.amber;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface AmberSample {
    byte DO_NOT_MATCH = (byte) 0;

    String sampleId();

    byte[] entries();

    default AmberPatient match(AmberSample other) {
        byte[] entries = entries();
        byte[] otherEntries = other.entries();

        if (entries.length != otherEntries.length) {
            throw new IllegalArgumentException("Unable to match different sized identities");
        }

        int match = 0;
        int validEntries = 0;
        for (int i = 0; i < entries.length; i++) {
            byte myByte = entries[i];
            byte otherByte = otherEntries[i];

            if (myByte != DO_NOT_MATCH && otherByte != DO_NOT_MATCH) {
                validEntries++;
                match += (myByte == otherByte ? 1 : 0);
            }
        }

        return ImmutableAmberPatient.builder()
                .sample(sampleId())
                .otherSample(other.sampleId())
                .sites(validEntries)
                .matches(match)
                .build();
    }

}
