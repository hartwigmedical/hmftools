package com.hartwig.hmftools.common.utils;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class ArrayDeckTest {

    @Test
    public void testGetterWorksAfterWraparound() {
        ArrayDeck<Integer> victim = new ArrayDeck<>(7);
        victim.add(0);
        for (int i = 1; i <= 5 ; i++) {
            victim.add(i);
            victim.removeFirst();
        }

        assertEquals(1, victim.size());
        assertEquals(5, victim.get(0).intValue());
        assertTrue(victim.tail >  victim.head);

        for (int i = 1; i <= 5; i++) {
            victim.add(5 + i);
            assertEquals(1 + i, victim.size());
            assertEquals(5 + i, victim.get(i).intValue());
        }

        assertTrue(victim.tail <  victim.head);
    }

}
