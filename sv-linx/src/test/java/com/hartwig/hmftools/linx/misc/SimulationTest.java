package com.hartwig.hmftools.linx.misc;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.linx.simulation.ShatteringConfig;
import com.hartwig.hmftools.linx.simulation.ShatteringResult;
import com.hartwig.hmftools.linx.simulation.ShatteringSim;

import org.junit.Test;

public class SimulationTest
{
    @Test
    public void testShatteringSim()
    {
        ShatteringConfig config = new ShatteringConfig(4, 1);
        ShatteringSim shatteringSim = new ShatteringSim(config, "");

        // shatteringSim.initialise(4, 1);

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

        ShatteringResult result = shatteringSim.getLatestResults();

        // assertEquals(0, results[SS_RESULTS_TEST_RUN]);
        assertEquals(4, result.segments());
        assertEquals(4, result.linkedSegments());
        assertEquals(5, result.exactMatches());
        assertEquals(5, result.adjacentSegments());

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

        result = shatteringSim.getLatestResults();

        // assertEquals(0, results[SS_RESULTS_TEST_RUN]);
        assertEquals(4, result.segments());
        assertEquals(4, result.linkedSegments());
        assertEquals(5, result.exactMatches());
        assertEquals(5, result.adjacentSegments());

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

        result = shatteringSim.getLatestResults();

        assertEquals(4, result.segments());
        assertEquals(2, result.linkedSegments());
        assertEquals(0, result.exactMatches());
        assertEquals(1, result.adjacentSegments());

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

        result = shatteringSim.getLatestResults();

        assertEquals(3, result.linkedSegments());
        assertEquals(3, result.exactMatches());
        assertEquals(3, result.adjacentSegments());
    }
}
