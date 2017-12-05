package com.hartwig.hmftools.common.purple.region;

import static com.hartwig.hmftools.common.purple.region.GermlineStatus.AMPLIFICATION;
import static com.hartwig.hmftools.common.purple.region.GermlineStatus.DIPLOID;
import static com.hartwig.hmftools.common.purple.region.GermlineStatus.HET_DELETION;
import static com.hartwig.hmftools.common.purple.region.GermlineStatus.HOM_DELETION;
import static com.hartwig.hmftools.common.purple.region.GermlineStatus.NOISE;
import static com.hartwig.hmftools.common.purple.region.GermlineStatus.UNKNOWN;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.purple.gender.Gender;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class GermlineStatusFactoryTest {

    private GermlineStatusFactory maleVictim = new GermlineStatusFactory(Gender.MALE);
    private GermlineStatusFactory femaleVictim = new GermlineStatusFactory(Gender.FEMALE);

    @Test
    public void testAutosome() {
        assertStatus("1", 0.00, UNKNOWN);
        assertStatus("1", 0.09, HOM_DELETION);
        assertStatus("1", 0.10, HET_DELETION);
        assertStatus("1", 0.79, HET_DELETION);
        assertStatus("1", 0.80, DIPLOID);
        assertStatus("1", 1.20, DIPLOID);
        assertStatus("1", 1.21, AMPLIFICATION);
        assertStatus("1", 2.20, AMPLIFICATION);
        assertStatus("1", 2.21, NOISE);
    }

    @Test
    public void testX() {
        assertStatus("X", 0.00, UNKNOWN);

        assertStatus("X", 0.049, HOM_DELETION, HOM_DELETION);
        assertStatus("X", 0.05, HET_DELETION, HOM_DELETION);

        assertStatus("X", 0.09, HET_DELETION, HOM_DELETION);
        assertStatus("X", 0.10, HET_DELETION, HET_DELETION);

        assertStatus("X", 0.39, HET_DELETION, HET_DELETION);
        assertStatus("X", 0.40, DIPLOID, HET_DELETION);

        assertStatus("X", 0.59, DIPLOID, HET_DELETION);
        assertStatus("X", 0.60, DIPLOID, HET_DELETION);
        assertStatus("X", 0.61, AMPLIFICATION, HET_DELETION);

        assertStatus("X", 0.79, AMPLIFICATION, HET_DELETION);
        assertStatus("X", 0.80, AMPLIFICATION, DIPLOID);

        assertStatus("X", 1.10, AMPLIFICATION, DIPLOID);
        assertStatus("X", 1.11, NOISE, DIPLOID);
        assertStatus("X", 1.21, NOISE, AMPLIFICATION);
        assertStatus("X", 2.20, NOISE, AMPLIFICATION);
        assertStatus("X", 2.21, NOISE, NOISE);
    }

    @Test
    public void testY() {
        assertStatus("Y", 0.00, UNKNOWN);
        assertStatus("Y", 0.049, HOM_DELETION);
        assertStatus("Y", 0.05, HET_DELETION);
        assertStatus("Y", 0.39, HET_DELETION);
        assertStatus("Y", 0.40, DIPLOID);
        assertStatus("Y", 0.60, DIPLOID);
        assertStatus("Y", 0.61, AMPLIFICATION);
        assertStatus("Y", 1.10, AMPLIFICATION);
        assertStatus("Y", 1.11, NOISE);
    }

    private void assertStatus(@NotNull final String chromosome, final double ratio, @NotNull final GermlineStatus expected) {
        assertEquals(expected, maleVictim.fromRatio(chromosome, ratio, 0.01));
        assertEquals(expected, femaleVictim.fromRatio(chromosome, ratio, 0.01));
    }

    private void assertStatus(@NotNull final String chromosome, final double ratio, @NotNull final GermlineStatus expectedMale,
            @NotNull final GermlineStatus expectedFemale) {
        assertEquals(expectedMale, maleVictim.fromRatio(chromosome, ratio, 0.01));
        assertEquals(expectedFemale, femaleVictim.fromRatio(chromosome, ratio, 0.01));
    }
}
