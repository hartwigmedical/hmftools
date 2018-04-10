package com.hartwig.hmftools.common.purple.qc;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.purple.gender.Gender;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class PurpleQCTest {

    @Test
    public void testSegmentScore() {
        assertEquals(PurpleQCStatus.PASS, create(240, 2d).status());
        assertEquals(PurpleQCStatus.FAIL_SEGMENT, create(241, 2d).status());

        assertEquals(PurpleQCStatus.PASS, create(120, 1d).status());
        assertEquals(PurpleQCStatus.FAIL_SEGMENT, create(121, 1d).status());
    }

    @Test
    public void testGenderCheck() {
        assertEquals(PurpleQCStatus.PASS, create(Gender.MALE, Gender.MALE).status());
        assertEquals(PurpleQCStatus.PASS, create(Gender.FEMALE, Gender.FEMALE).status());
        assertEquals(PurpleQCStatus.PASS, create(Gender.FEMALE, Gender.MALE_KLINEFELTER).status());

        assertEquals(PurpleQCStatus.FAIL_GENDER, create(Gender.MALE, Gender.FEMALE).status());
        assertEquals(PurpleQCStatus.FAIL_GENDER, create(Gender.FEMALE, Gender.MALE).status());
        assertEquals(PurpleQCStatus.FAIL_GENDER, create(Gender.MALE, Gender.MALE_KLINEFELTER).status());
    }

    @NotNull
    private static PurpleQC create(int unsupportedSegments, double ploidy) {
        return ImmutablePurpleQC.builder()
                .unsupportedSegments(unsupportedSegments)
                .ploidy(ploidy)
                .amberGender(Gender.MALE)
                .cobaltGender(Gender.MALE)
                .build();
    }

    @NotNull
    private static PurpleQC create(@NotNull final Gender amberGender, @NotNull final Gender cobaltGender) {
        return ImmutablePurpleQC.builder().unsupportedSegments(1).ploidy(1).amberGender(amberGender).cobaltGender(cobaltGender).build();
    }
}
