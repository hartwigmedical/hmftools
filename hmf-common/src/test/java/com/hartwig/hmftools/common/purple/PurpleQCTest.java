package com.hartwig.hmftools.common.purple;

import static com.hartwig.hmftools.common.purple.PurpleQCStatus.MAX_DELETED_GENES;
import static com.hartwig.hmftools.common.purple.PurpleQCStatus.PASS;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.chromosome.GermlineAberration;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class PurpleQCTest
{
    @Test
    public void testDefault()
    {
        ImmutablePurpleQC.Builder template = builder();
        // assertTrue(template.build().pass());
        assertStatus(template, PASS);
    }

    @Test
    public void testGender()
    {
        ImmutablePurpleQC.Builder template = builder().amberGender(Gender.FEMALE).cobaltGender(Gender.FEMALE);
        assertStatus(template, PASS);
        assertStatus(template.cobaltGender(Gender.MALE), PurpleQCStatus.WARN_GENDER_MISMATCH);
        assertStatus(template.addGermlineAberrations(GermlineAberration.KLINEFELTER), PASS);
    }

    @Test
    public void testNoTumor()
    {
        assertStatus(builder().method(FittedPurityMethod.NO_TUMOR), PurpleQCStatus.FAIL_NO_TUMOR);
    }

    @Test
    public void testIgnoreLowPurityIfNoTumor()
    {
        assertStatus(builder().purity(0.01).method(FittedPurityMethod.NO_TUMOR), PurpleQCStatus.FAIL_NO_TUMOR);
    }

    @Test
    public void testContamination()
    {
        assertStatus(builder().contamination(PurpleQCStatus.MAX_CONTAMINATION), PASS);
        assertStatus(builder().contamination(PurpleQCStatus.MAX_CONTAMINATION + 0.01), PurpleQCStatus.FAIL_CONTAMINATION);
    }

    @Test
    public void testDeletedGenes()
    {
        assertStatus(builder().deletedGenes(PurpleQCStatus.MAX_DELETED_GENES), PASS);
        assertStatus(builder().deletedGenes(PurpleQCStatus.MAX_DELETED_GENES + 1), PurpleQCStatus.WARN_DELETED_GENES);
    }

    @Test
    public void testHighCopyNumberNoise()
    {
        assertStatus(builder().unsupportedCopyNumberSegments(PurpleQCStatus.MAX_UNSUPPORTED_SEGMENTS), PASS);
        assertStatus(builder().unsupportedCopyNumberSegments(PurpleQCStatus.MAX_UNSUPPORTED_SEGMENTS + 1),
                PurpleQCStatus.WARN_HIGH_COPY_NUMBER_NOISE);
    }

    @Test
    public void testLowPurity()
    {
        assertStatus(builder().purity(PurpleQCStatus.MIN_PURITY), PASS);
        assertStatus(builder().purity(PurpleQCStatus.MIN_PURITY - 0.01), PurpleQCStatus.WARN_LOW_PURITY);
    }

    @Test
    public void testMultiple()
    {
        ImmutablePurpleQC.Builder victim = builder()
            .purity(PurpleQCStatus.MIN_PURITY - 0.01)
            .contamination(PurpleQCStatus.MAX_CONTAMINATION + 0.01);

        Set<PurpleQCStatus> requiredStatusSet = Sets.newHashSet(PurpleQCStatus.WARN_LOW_PURITY);
        requiredStatusSet.add(PurpleQCStatus.FAIL_CONTAMINATION);
        assertStatus(victim, requiredStatusSet);
    }

    private static void assertStatus(final ImmutablePurpleQC.Builder victim, final PurpleQCStatus status)
    {
        Set<PurpleQCStatus> requiredStatusSet = Sets.newHashSet(status);
        assertStatus(victim, requiredStatusSet);
    }

    private static void assertStatus(final ImmutablePurpleQC.Builder victim, final Set<PurpleQCStatus> requiredStatusSet)
    {
        PurpleQC qc = victim.status(Sets.newHashSet(PASS)).build();

        Set<PurpleQCStatus> statusSet = PurpleQCStatus.calcStatus(
                PurpleQCStatus.genderPass(qc.amberGender(), qc.cobaltGender(), qc.germlineAberrations()),
                qc.unsupportedCopyNumberSegments(), qc.deletedGenes(), qc.purity(), qc.method(), qc.contamination(), MAX_DELETED_GENES);

        assertTrue(statusSet.size() == requiredStatusSet.size());

        assertTrue(requiredStatusSet.stream().allMatch(x -> statusSet.contains(x)));
    }

    /*
    private static void assertStatus(@NotNull PurpleQC victim, @NotNull PurpleQCStatus... status)
    {
        assertEquals(victim.status().size(), status.length);
        for(PurpleQCStatus purpleQCStatus : status)
        {
            assertTrue(victim.status().contains(purpleQCStatus));
        }
    }
    */

    @NotNull
    private static ImmutablePurpleQC.Builder builder()
    {
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
