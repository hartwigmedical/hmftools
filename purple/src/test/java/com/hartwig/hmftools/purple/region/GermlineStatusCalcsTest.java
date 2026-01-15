package com.hartwig.hmftools.purple.region;

import static com.hartwig.hmftools.common.purple.GermlineStatus.AMPLIFICATION;
import static com.hartwig.hmftools.common.purple.GermlineStatus.DIPLOID;
import static com.hartwig.hmftools.common.purple.GermlineStatus.HET_DELETION;
import static com.hartwig.hmftools.common.purple.GermlineStatus.HOM_DELETION;
import static com.hartwig.hmftools.common.purple.GermlineStatus.LIKELY_DIPLOID;
import static com.hartwig.hmftools.common.purple.GermlineStatus.NOISE;
import static com.hartwig.hmftools.common.purple.GermlineStatus.UNKNOWN;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.cobalt.CobaltTestUtils;
import com.hartwig.hmftools.common.purple.GermlineStatus;

import org.junit.Test;

public class GermlineStatusCalcsTest
{
    private final GermlineStatusCalcs mMaleCalcs = new GermlineStatusCalcs(CobaltTestUtils.male());
    private final GermlineStatusCalcs mFemaleCalcs = new GermlineStatusCalcs(CobaltTestUtils.female());

    @Test
    public void testAutosome()
    {
        assertStatus("1", 0.00, HOM_DELETION);
        assertStatus("1", 0.09, HOM_DELETION);
        assertStatus("1", 0.10, HET_DELETION);
        assertStatus("1", 0.69, HET_DELETION);
        assertStatus("1", 0.70, LIKELY_DIPLOID);
        assertStatus("1", 0.84, LIKELY_DIPLOID);
        assertStatus("1", 0.8499, LIKELY_DIPLOID);
        assertStatus("1", 0.85, DIPLOID);
        assertStatus("1", 0.85001, DIPLOID);
        assertStatus("1", 1.1499, DIPLOID);
        assertStatus("1", 1.15, LIKELY_DIPLOID);
        assertStatus("1", 1.15001, LIKELY_DIPLOID);
        assertStatus("1", 1.16, LIKELY_DIPLOID);
        assertStatus("1", 1.299, LIKELY_DIPLOID);
        assertStatus("1", 1.3, AMPLIFICATION);
        assertStatus("1", 1.3001, AMPLIFICATION);
        assertStatus("1", 2.199, AMPLIFICATION);
        assertStatus("1", 2.20, NOISE);
        assertStatus("1", 2.21, NOISE);
    }

    @Test
    public void testX()
    {
        assertStatus("X", 0.00, HOM_DELETION);

        assertStatus("X", 0.049, HOM_DELETION, HOM_DELETION);
        assertStatus("X", 0.05, HET_DELETION, HOM_DELETION);

        assertStatus("X", 0.09, HET_DELETION, HOM_DELETION);
        assertStatus("X", 0.10, HET_DELETION, HET_DELETION);

        assertStatus("X", 0.424, LIKELY_DIPLOID, HET_DELETION);
        assertStatus("X", 0.425, DIPLOID, HET_DELETION);

        assertStatus("X", 0.574, DIPLOID, HET_DELETION);
        assertStatus("X", 0.575, LIKELY_DIPLOID, HET_DELETION);
        assertStatus("X", 0.576, LIKELY_DIPLOID, HET_DELETION);

        assertStatus("X", 0.84, AMPLIFICATION, LIKELY_DIPLOID);
        assertStatus("X", 0.85, AMPLIFICATION, DIPLOID);

        assertStatus("X", 1.10, NOISE, DIPLOID);
        assertStatus("X", 1.15, NOISE, LIKELY_DIPLOID);
        assertStatus("X", 1.16, NOISE, LIKELY_DIPLOID);
        assertStatus("X", 1.8, NOISE, AMPLIFICATION);
        assertStatus("X", 2.20, NOISE, NOISE);
        assertStatus("X", 2.21, NOISE, NOISE);
    }

    @Test
    public void testY()
    {
        assertStatus("Y", 0.00, HOM_DELETION);
        assertStatus("Y", 0.049, HOM_DELETION);
        assertStatus("Y", 0.05, HET_DELETION);
        assertStatus("Y", 0.424, LIKELY_DIPLOID);
        assertStatus("Y", 0.425, DIPLOID);
        assertStatus("Y", 0.575, LIKELY_DIPLOID);
        assertStatus("Y", 0.576, LIKELY_DIPLOID);
        assertStatus("Y", 0.9, AMPLIFICATION);
        assertStatus("Y", 1.10, NOISE);
        assertStatus("Y", 1.11, NOISE);
    }

    private void assertStatus(final String chromosome, final double ratio, final GermlineStatus expected)
    {
        assertEquals(expected, mMaleCalcs.calcStatus(chromosome, ratio, 0.01, 1));
        if(chromosome.equals("Y"))
        {
            assertEquals(UNKNOWN, mFemaleCalcs.calcStatus(chromosome, ratio, 0.01, 1));
        }
        else
        {
            assertEquals(expected, mFemaleCalcs.calcStatus(chromosome, ratio, 0.01, 1));
        }
    }

    private void assertStatus(
            final String chromosome, final double ratio, final GermlineStatus expectedMale, final GermlineStatus expectedFemale)
    {
        assertEquals(expectedMale, mMaleCalcs.calcStatus(chromosome, ratio, 0.01, 1));
        assertEquals(expectedFemale, mFemaleCalcs.calcStatus(chromosome, ratio, 0.01, 1));
    }
}
