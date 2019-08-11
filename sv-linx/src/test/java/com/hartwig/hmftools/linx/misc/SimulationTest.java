package com.hartwig.hmftools.linx.misc;

import static com.hartwig.hmftools.linx.simulation.SvSimShattering.SS_RESULTS_ADJACENT_SEGS;
import static com.hartwig.hmftools.linx.simulation.SvSimShattering.SS_RESULTS_EXACT_MATCHES;
import static com.hartwig.hmftools.linx.simulation.SvSimShattering.SS_RESULTS_SEGMENTS;
import static com.hartwig.hmftools.linx.simulation.SvSimShattering.SS_RESULTS_SEGS_LINKED;
import static com.hartwig.hmftools.linx.simulation.SvSimShattering.SS_RESULTS_TEST_RUN;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.linx.simulation.SvSimShattering;

import org.junit.Test;

public class SimulationTest
{
    @Test
    public void testShatteringSim()
    {
        SvSimShattering shatteringSim = new SvSimShattering();

        shatteringSim.initialise(4, 1);

        // first test all segments added first to last

        // linking order
        List<Integer> linkOrder = Lists.newArrayList();
        linkOrder.add(0);
        linkOrder.add(0);
        linkOrder.add(0);
        linkOrder.add(0);
        linkOrder.add(0);
        linkOrder.add(0);
        linkOrder.add(0);
        linkOrder.add(0);
        linkOrder.add(0);
        linkOrder.add(0);

        shatteringSim.setSpecifiedOrder(linkOrder);

        shatteringSim.run();
        assertTrue(shatteringSim.validRun());

        int[] results = shatteringSim.getLatestResults();

        // assertEquals(0, results[SS_RESULTS_TEST_RUN]);
        assertEquals(4, results[SS_RESULTS_SEGMENTS]);
        assertEquals(4, results[SS_RESULTS_SEGS_LINKED]);
        assertEquals(5, results[SS_RESULTS_EXACT_MATCHES]);
        assertEquals(5, results[SS_RESULTS_ADJACENT_SEGS]);

        // test all segments added last to first

        linkOrder.clear();

        linkOrder.add(9);
        linkOrder.add(8);
        linkOrder.add(7);
        linkOrder.add(6);
        linkOrder.add(5);
        linkOrder.add(4);
        linkOrder.add(3);
        linkOrder.add(2);
        linkOrder.add(1);
        linkOrder.add(0);

        shatteringSim.setSpecifiedOrder(linkOrder);

        shatteringSim.run();
        assertTrue(shatteringSim.validRun());

        results = shatteringSim.getLatestResults();

        // assertEquals(0, results[SS_RESULTS_TEST_RUN]);
        assertEquals(4, results[SS_RESULTS_SEGMENTS]);
        assertEquals(4, results[SS_RESULTS_SEGS_LINKED]);
        assertEquals(5, results[SS_RESULTS_EXACT_MATCHES]);
        assertEquals(5, results[SS_RESULTS_ADJACENT_SEGS]);

        // test only 2 segs added
        linkOrder.clear();

        linkOrder.add(0);
        linkOrder.add(3);
        linkOrder.add(3);
        linkOrder.add(7);
        linkOrder.add(2);
        linkOrder.add(3);

        shatteringSim.setSpecifiedOrder(linkOrder);

        shatteringSim.run();
        assertTrue(shatteringSim.validRun());

        results = shatteringSim.getLatestResults();

        assertEquals(4, results[SS_RESULTS_SEGMENTS]);
        assertEquals(2, results[SS_RESULTS_SEGS_LINKED]);
        assertEquals(0, results[SS_RESULTS_EXACT_MATCHES]);
        assertEquals(1, results[SS_RESULTS_ADJACENT_SEGS]);

        // leave a single fully open link
        linkOrder.clear();

        linkOrder.add(9);
        linkOrder.add(8);
        linkOrder.add(7);
        linkOrder.add(6);
        linkOrder.add(0);
        linkOrder.add(1);
        linkOrder.add(0);
        linkOrder.add(3);

        shatteringSim.setSpecifiedOrder(linkOrder);

        shatteringSim.run();
        assertTrue(shatteringSim.validRun());

        results = shatteringSim.getLatestResults();

        // assertEquals(0, results[SS_RESULTS_TEST_RUN]);
        assertEquals(3, results[SS_RESULTS_SEGS_LINKED]);
        assertEquals(3, results[SS_RESULTS_EXACT_MATCHES]);
        assertEquals(3, results[SS_RESULTS_ADJACENT_SEGS]);


    }
}
