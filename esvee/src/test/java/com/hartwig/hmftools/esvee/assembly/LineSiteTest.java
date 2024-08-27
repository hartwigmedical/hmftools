package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.esvee.common.CommonUtils.withLineProximity;
import static com.hartwig.hmftools.esvee.common.SvConstants.LINE_INDEL_MAX_GAP;
import static com.hartwig.hmftools.esvee.common.SvConstants.LINE_INDEL_MAX_OVERLAP;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class LineSiteTest
{
    @Test
    public void testLineInsertionProximity()
    {
        assertFalse(withLineProximity(100, 100, FORWARD, FORWARD));

        assertTrue(withLineProximity(100, 100, FORWARD, REVERSE));
        assertTrue(withLineProximity(100, 100, REVERSE, FORWARD));

        // DEL within range
        int lowerPosition = 100;
        int upperPosition = lowerPosition + LINE_INDEL_MAX_GAP - 1;
        assertTrue(withLineProximity(lowerPosition, upperPosition, FORWARD, REVERSE));
        assertTrue(withLineProximity(upperPosition, lowerPosition, REVERSE, FORWARD));

        // too far for a DEL
        assertTrue(withLineProximity(lowerPosition, upperPosition + 1, FORWARD, REVERSE));
        assertTrue(withLineProximity(upperPosition + 1, lowerPosition, REVERSE, FORWARD));

        // DUP within range
        upperPosition = lowerPosition + LINE_INDEL_MAX_OVERLAP - 1;
        assertTrue(withLineProximity(lowerPosition, upperPosition, REVERSE, FORWARD));
        assertTrue(withLineProximity(upperPosition, lowerPosition, FORWARD, REVERSE));

        // too far for a DUP
        assertTrue(withLineProximity(lowerPosition, upperPosition + 1, REVERSE, FORWARD));
        assertTrue(withLineProximity(upperPosition + 1, lowerPosition, FORWARD, REVERSE));
    }

}
