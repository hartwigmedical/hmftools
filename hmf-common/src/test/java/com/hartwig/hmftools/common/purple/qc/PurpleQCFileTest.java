package com.hartwig.hmftools.common.purple.qc;

import static org.junit.Assert.assertEquals;

import java.io.IOException;

import com.hartwig.hmftools.common.purple.gender.Gender;

import org.junit.Test;

public class PurpleQCFileTest {

    @Test
    public void testReadWrite() throws IOException {
        final PurpleQC expected = create();
        assertEquals(expected, PurpleQCFile.fromLines(PurpleQCFile.toLines(expected)));
    }

    private static PurpleQC create() {
        return ImmutablePurpleQC.builder()
                .trailingSegments(90)
                .ratioSegments(34)
                .ploidy(3.17)
                .purpleGender(Gender.MALE)
                .cobaltGender(Gender.FEMALE)
                .build();
    }

}
