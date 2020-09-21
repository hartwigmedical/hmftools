package com.hartwig.hmftools.common.purple.qc;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.genome.chromosome.GermlineAberration;
import com.hartwig.hmftools.common.purple.gender.Gender;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class PurpleQCTest {

    @Test
    public void testSegmentScore() {
        assertEquals(PurpleQCStatus.PASS, create(220, 2d).status());
        assertEquals(PurpleQCStatus.FAIL_SEGMENT, create(221, 2d).status());

        assertEquals(PurpleQCStatus.PASS, create(220, 1d).status());
        assertEquals(PurpleQCStatus.FAIL_SEGMENT, create(221, 1d).status());
    }

    @Test
    public void testGenderCheck() {
        assertEquals(PurpleQCStatus.PASS, create(Gender.MALE, Gender.MALE).status());
        assertEquals(PurpleQCStatus.PASS, create(Gender.FEMALE, Gender.FEMALE).status());

        assertEquals(PurpleQCStatus.FAIL_GENDER, create(Gender.MALE, Gender.FEMALE).status());
        assertEquals(PurpleQCStatus.FAIL_GENDER, create(Gender.FEMALE, Gender.MALE).status());

        assertEquals(PurpleQCStatus.PASS, create(Gender.FEMALE, Gender.MALE, GermlineAberration.KLINEFELTER).status());
        assertEquals(PurpleQCStatus.PASS, create(Gender.MALE, Gender.MALE).status());
    }

    @Test
    public void testDeletedGenes() {
        assertEquals(PurpleQCStatus.PASS, create(280).status());
        assertEquals(PurpleQCStatus.FAIL_DELETED_GENES, create(281).status());
    }

    @NotNull
    private static PurpleQC create(int deletedGenes) {
        return ImmutablePurpleQC.builder()
                .unsupportedSegments(1)
                .ploidy(2)
                .amberGender(Gender.MALE)
                .cobaltGender(Gender.MALE)
                .deletedGenes(deletedGenes)
                .build();
    }

    @NotNull
    private static PurpleQC create(int unsupportedSegments, double ploidy) {
        return ImmutablePurpleQC.builder()
                .unsupportedSegments(unsupportedSegments)
                .ploidy(ploidy)
                .amberGender(Gender.MALE)
                .cobaltGender(Gender.MALE)
                .deletedGenes(1)
                .build();
    }

    @NotNull
    private static PurpleQC create(@NotNull final Gender amberGender, @NotNull final Gender cobaltGender,
            @NotNull final GermlineAberration... aberrations) {
        return ImmutablePurpleQC.builder()
                .unsupportedSegments(1)
                .ploidy(1)
                .amberGender(amberGender)
                .cobaltGender(cobaltGender)
                .deletedGenes(1)
                .addGermlineAberrations(aberrations)
                .build();
    }
}
