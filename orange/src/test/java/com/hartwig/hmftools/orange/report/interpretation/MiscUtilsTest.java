package com.hartwig.hmftools.orange.report.interpretation;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.orange.report.tables.TableCommon;

import org.junit.Test;

public class MiscUtilsTest
{
    @Test
    public void canPrefixChromosomes()
    {
        assertEquals("05", TableCommon.zeroPrefixed("5"));
        assertEquals("05p13.2", TableCommon.zeroPrefixed("5p13.2"));

        assertEquals("15", TableCommon.zeroPrefixed("15"));
        assertEquals("15q21.1", TableCommon.zeroPrefixed("15q21.1"));

        assertEquals("X", TableCommon.zeroPrefixed("X"));
    }
}