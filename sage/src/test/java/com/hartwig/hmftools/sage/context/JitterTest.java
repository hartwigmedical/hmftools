package com.hartwig.hmftools.sage.context;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class JitterTest {

    @Test
    public void testPolyA() {
        String shorter = "GATCAAAAAAAAAGATC";
        String ref = "GATCAAAAAAAAAAGATC";
        String longer = "GATCAAAAAAAAAAAGATC";

        RealignedContext result = null;

        result = Realigned.jitter(3, 2, 16, ref.getBytes(), 3, ref.getBytes());
        assertRealigned(result, RealignedType.EXACT, 0);

        result = Realigned.jitter(3, 2, 16, ref.getBytes(), 3, shorter.getBytes());
        assertRealigned(result, RealignedType.SHORTENED, 9);

        result = Realigned.jitter(3, 2, 16, ref.getBytes(), 3, longer.getBytes());
        assertRealigned(result, RealignedType.LENGTHENED, 11);
    }

    @Test
    public void testDiNucleotideRepeat() {
        String shorter = "GATCATATATATGATC";
        String ref = "GATCATATATATATGATC";
        String longer = "GATCATATATATATATGATC";

        RealignedContext result = null;
        result = Realigned.jitter(3, 2, 16, ref.getBytes(), 3, ref.getBytes());
        assertRealigned(result, RealignedType.EXACT, 0);

        result = Realigned.jitter(3, 2, 16, ref.getBytes(), 3, shorter.getBytes());
        assertRealigned(result, RealignedType.SHORTENED, 4);

        result = Realigned.jitter(3, 2, 16, ref.getBytes(), 3, longer.getBytes());
        assertRealigned(result, RealignedType.LENGTHENED, 6);
    }

    private void assertRealigned(RealignedContext victim, RealignedType expectedType, int expectedRepeatCount) {
        assertEquals(expectedType, victim.type());
        assertEquals(expectedRepeatCount, victim.repeatCount());
    }

}
