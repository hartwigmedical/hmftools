package com.hartwig.hmftools.common.purple.region;

import static com.hartwig.hmftools.common.purple.region.ObservedRegionStatus.GERMLINE_AMPLIFICATION;
import static com.hartwig.hmftools.common.purple.region.ObservedRegionStatus.GERMLINE_HET_DELETION;
import static com.hartwig.hmftools.common.purple.region.ObservedRegionStatus.GERMLINE_HOM_DELETION;
import static com.hartwig.hmftools.common.purple.region.ObservedRegionStatus.GERMLINE_NOISE;
import static com.hartwig.hmftools.common.purple.region.ObservedRegionStatus.SOMATIC;
import static com.hartwig.hmftools.common.purple.region.ObservedRegionStatus.UNKNOWN;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.purple.gender.Gender;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ObservedRegionStatusFactoryTest {

    private ObservedRegionStatusFactory maleVictim = new ObservedRegionStatusFactory(Gender.MALE);
    private ObservedRegionStatusFactory femaleVictim = new ObservedRegionStatusFactory(Gender.FEMALE);

    @Test
    public void testAutosome() {
        assertStatus("1", 0.00, UNKNOWN);
        assertStatus("1", 0.09, GERMLINE_HOM_DELETION);
        assertStatus("1", 0.10, GERMLINE_HET_DELETION);
        assertStatus("1", 0.79, GERMLINE_HET_DELETION);
        assertStatus("1", 0.80, SOMATIC);
        assertStatus("1", 1.20, SOMATIC);
        assertStatus("1", 1.21, GERMLINE_AMPLIFICATION);
        assertStatus("1", 2.20, GERMLINE_AMPLIFICATION);
        assertStatus("1", 2.21, GERMLINE_NOISE);
    }

    @Test
    public void testX() {
        assertStatus("X", 0.00, UNKNOWN);

        assertStatus("X", 0.049, GERMLINE_HOM_DELETION, GERMLINE_HOM_DELETION);
        assertStatus("X", 0.05, GERMLINE_HET_DELETION, GERMLINE_HOM_DELETION);

        assertStatus("X", 0.09, GERMLINE_HET_DELETION, GERMLINE_HOM_DELETION);
        assertStatus("X", 0.10, GERMLINE_HET_DELETION, GERMLINE_HET_DELETION);

        assertStatus("X", 0.39, GERMLINE_HET_DELETION, GERMLINE_HET_DELETION);
        assertStatus("X", 0.40, SOMATIC, GERMLINE_HET_DELETION);

        assertStatus("X", 0.59, SOMATIC, GERMLINE_HET_DELETION);
        assertStatus("X", 0.60, SOMATIC, GERMLINE_HET_DELETION);
        assertStatus("X", 0.61, GERMLINE_AMPLIFICATION, GERMLINE_HET_DELETION);

        assertStatus("X", 0.79, GERMLINE_AMPLIFICATION, GERMLINE_HET_DELETION);
        assertStatus("X", 0.80, GERMLINE_AMPLIFICATION, SOMATIC);

        assertStatus("X", 1.10, GERMLINE_AMPLIFICATION, SOMATIC);
        assertStatus("X", 1.11, GERMLINE_NOISE, SOMATIC);
        assertStatus("X", 1.21, GERMLINE_NOISE, GERMLINE_AMPLIFICATION);
        assertStatus("X", 2.20, GERMLINE_NOISE, GERMLINE_AMPLIFICATION);
        assertStatus("X", 2.21, GERMLINE_NOISE, GERMLINE_NOISE);
    }

    @Test
    public void testY() {
        assertStatus("Y", 0.00, UNKNOWN);
        assertStatus("Y", 0.049, GERMLINE_HOM_DELETION);
        assertStatus("Y", 0.05, GERMLINE_HET_DELETION);
        assertStatus("Y", 0.39, GERMLINE_HET_DELETION);
        assertStatus("Y", 0.40, SOMATIC);
        assertStatus("Y", 0.60, SOMATIC);
        assertStatus("Y", 0.61, GERMLINE_AMPLIFICATION);
        assertStatus("Y", 1.10, GERMLINE_AMPLIFICATION);
        assertStatus("Y", 1.11, GERMLINE_NOISE);
    }

    private void assertStatus(@NotNull final String chromosome, final double ratio, @NotNull final ObservedRegionStatus expected) {
        assertEquals(expected, maleVictim.fromRatio(chromosome, ratio, 0.01));
        assertEquals(expected, femaleVictim.fromRatio(chromosome, ratio, 0.01));

    }

    private void assertStatus(@NotNull final String chromosome, final double ratio, @NotNull final ObservedRegionStatus expectedMale,
            @NotNull final ObservedRegionStatus expectedFemale) {
        assertEquals(expectedMale, maleVictim.fromRatio(chromosome, ratio, 0.01));
        assertEquals(expectedFemale, femaleVictim.fromRatio(chromosome, ratio, 0.01));

    }

}
