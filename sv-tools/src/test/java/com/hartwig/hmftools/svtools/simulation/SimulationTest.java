package com.hartwig.hmftools.svtools.simulation;

import static com.hartwig.hmftools.svtools.simulation.ShatteringSim.calcLinkCount;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;

import org.junit.Test;

public class SimulationTest
{
    @Test
    public void testShatteringSim()
    {
        int segCount = 4;
        ShatteringConfig config = new ShatteringConfig(segCount, 1);
        ShatteringSim shatteringSim = new ShatteringSim(config, "");

        // Configurator.setRootLevel(Level.DEBUG);

        // first test all segments added first to last

        // linking order
        List<Integer> linkOrder = Lists.newArrayList();
        int linkCount = calcLinkCount(segCount + 2);

        for(int i = 0; i <= linkCount; ++i)
        {
            linkOrder.add(i);
        }

        shatteringSim.setSpecifiedOrder(linkOrder);

        shatteringSim.run();
        assertTrue(shatteringSim.validRun());

        ShatteringResult result = shatteringSim.getLatestResults();

        // assertEquals(0, results[SS_RESULTS_TEST_RUN]);
        assertEquals(4, result.segments());
        assertEquals(4, result.linkedSegments());
        assertEquals(5, result.exactRepairs());
        assertEquals(5, result.adjacentSegments());
        assertEquals(0, result.inferredLinks());
        assertEquals(0, result.inferredLost());

        // test all segments added last to first
        linkOrder.clear();

        for(int i = linkCount; i >= 0; --i)
        {
            linkOrder.add(i);
        }

        shatteringSim.setSpecifiedOrder(linkOrder);

        shatteringSim.run();
        assertTrue(shatteringSim.validRun());

        result = shatteringSim.getLatestResults();

        assertEquals(4, result.segments());
        assertEquals(4, result.linkedSegments());
        assertEquals(5, result.exactRepairs());
        assertEquals(5, result.adjacentSegments());
        assertEquals(0, result.inferredLinks());
        assertEquals(0, result.inferredLost());

        // test no segments added
        linkOrder.clear();
        linkOrder.add(0);
        linkOrder.add(9);

        shatteringSim.setSpecifiedOrder(linkOrder);

        shatteringSim.run();
        assertTrue(shatteringSim.validRun());

        result = shatteringSim.getLatestResults();

        assertEquals(0, result.linkedSegments());
        assertEquals(0, result.exactRepairs());
        assertEquals(0, result.adjacentSegments());
        assertEquals(1, result.inferredLost());

        // test only 2 segs added
        linkOrder.clear();

        linkOrder.add(0);
        linkOrder.add(3);
        linkOrder.add(4);
        linkOrder.add(7);
        linkOrder.add(8);
        linkOrder.add(9);

        shatteringSim.setSpecifiedOrder(linkOrder);

        shatteringSim.run();
        assertTrue(shatteringSim.validRun());

        result = shatteringSim.getLatestResults();

        assertEquals(4, result.segments());
        assertEquals(2, result.linkedSegments());
        assertEquals(1, result.exactRepairs());
        assertEquals(1, result.adjacentSegments());
        assertEquals(1, result.inferredLinks());
        assertEquals(2, result.inferredLost());

        // leave a single fully open link
        linkOrder.clear();

        linkOrder.add(9);
        linkOrder.add(7);
        linkOrder.add(8); // reversed
        linkOrder.add(6);
        linkOrder.add(5);
        linkOrder.add(2);
        linkOrder.add(1);
        linkOrder.add(0);

        shatteringSim.setSpecifiedOrder(linkOrder);

        shatteringSim.run();
        assertTrue(shatteringSim.validRun());

        result = shatteringSim.getLatestResults();

        assertEquals(3, result.linkedSegments());
        assertEquals(1, result.exactRepairs());
        assertEquals(3, result.adjacentSegments());
        assertEquals(2, result.inferredLinks());
        assertEquals(1, result.inferredLost());
    }
}
