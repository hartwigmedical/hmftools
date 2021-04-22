package com.hartwig.hmftools.common.purple.qc;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.genome.chromosome.GermlineAberration;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.purity.FittedPurityMethod;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class PurpleQCTest {

    @Test
    public void testDefault() {
        ImmutablePurpleQC.Builder template = builder();
        assertTrue(template.build().pass());
        assertStatus(template, PurpleQCStatus.PASS);
    }

    @Test
    public void testOldGender() {
        ImmutablePurpleQC.Builder template = builder().amberGender(Gender.FEMALE).cobaltGender(Gender.MALE_KLINEFELTER);
        assertStatus(template, PurpleQCStatus.PASS);
    }

    @Test
    public void testGender() {
        ImmutablePurpleQC.Builder template = builder().amberGender(Gender.FEMALE).cobaltGender(Gender.FEMALE);
        assertStatus(template, PurpleQCStatus.PASS);
        assertStatus(template.cobaltGender(Gender.MALE), PurpleQCStatus.WARN_GENDER_MISMATCH);
        assertStatus(template.addGermlineAberrations(GermlineAberration.KLINEFELTER), PurpleQCStatus.PASS);
    }

    @Test
    public void testNoTumor() {
        assertStatus(builder().method(FittedPurityMethod.NO_TUMOR), PurpleQCStatus.FAIL_NO_TUMOR);
    }

    @Test
    public void testIgnoreLowPurityIfNoTumor() {
        assertStatus(builder().purity(0.01).method(FittedPurityMethod.NO_TUMOR), PurpleQCStatus.FAIL_NO_TUMOR);
    }

    @Test
    public void testContamination() {
        assertStatus(builder().contamination(PurpleQC.MAX_CONTAMINATION), PurpleQCStatus.PASS);
        assertStatus(builder().contamination(PurpleQC.MAX_CONTAMINATION + 0.01), PurpleQCStatus.FAIL_CONTAMINATION);
    }

    @Test
    public void testDeletedGenes() {
        assertStatus(builder().deletedGenes(PurpleQC.MAX_DELETED_GENES), PurpleQCStatus.PASS);
        assertStatus(builder().deletedGenes(PurpleQC.MAX_DELETED_GENES + 1), PurpleQCStatus.WARN_DELETED_GENES);
    }

    @Test
    public void testHighCopyNumberNoise() {
        assertStatus(builder().unsupportedCopyNumberSegments(PurpleQC.MAX_UNSUPPORTED_SEGMENTS), PurpleQCStatus.PASS);
        assertStatus(builder().unsupportedCopyNumberSegments(PurpleQC.MAX_UNSUPPORTED_SEGMENTS + 1),
                PurpleQCStatus.WARN_HIGH_COPY_NUMBER_NOISE);
    }

    @Test
    public void testLowPurity() {
        assertStatus(builder().purity(PurpleQC.MIN_PURITY), PurpleQCStatus.PASS);
        assertStatus(builder().purity(PurpleQC.MIN_PURITY - 0.01), PurpleQCStatus.WARN_LOW_PURITY);
    }

    @Test
    public void testMultiple() {
        PurpleQC victim = builder().purity(PurpleQC.MIN_PURITY - 0.01).contamination(PurpleQC.MAX_CONTAMINATION + 0.01).build();
        assertStatus(victim, PurpleQCStatus.WARN_LOW_PURITY, PurpleQCStatus.FAIL_CONTAMINATION);
    }

    private static void assertStatus(@NotNull ImmutablePurpleQC.Builder victim, @NotNull PurpleQCStatus... status) {
        assertStatus(victim.build(), status);
    }

    private static void assertStatus(@NotNull PurpleQC victim, @NotNull PurpleQCStatus... status) {
        assertEquals(victim.status().size(), status.length);
        for (PurpleQCStatus purpleQCStatus : status) {
            assertTrue(victim.status().contains(purpleQCStatus));
        }
    }

    @NotNull
    private static ImmutablePurpleQC.Builder builder() {
        return ImmutablePurpleQC.builder()
                .contamination(0)
                .copyNumberSegments(100)
                .unsupportedCopyNumberSegments(0)
                .purity(1)
                .amberGender(Gender.MALE)
                .cobaltGender(Gender.MALE)
                .deletedGenes(0)
                .method(FittedPurityMethod.NORMAL)
                .amberMeanDepth(91);
    }
}
