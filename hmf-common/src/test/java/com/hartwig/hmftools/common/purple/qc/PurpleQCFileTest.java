package com.hartwig.hmftools.common.purple.qc;

import static org.junit.Assert.assertEquals;

import java.io.IOException;

import com.hartwig.hmftools.common.purple.gender.Gender;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class PurpleQCFileTest {

    @Test
    public void testReadWrite() throws IOException {
        final PurpleQC expected = create();
        assertEquals(expected, PurpleQCFile.fromLines(PurpleQCFile.toLines(expected)));
    }

    @NotNull
    private static PurpleQC create() {
        return ImmutablePurpleQC.builder()
                .unsupportedSegments(90)
                .ploidy(3.17)
                .amberGender(Gender.MALE)
                .cobaltGender(Gender.FEMALE)
                .build();
    }
}
