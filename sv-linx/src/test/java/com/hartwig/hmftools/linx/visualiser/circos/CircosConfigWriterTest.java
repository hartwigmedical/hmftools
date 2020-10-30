package com.hartwig.hmftools.linx.visualiser.circos;

import static org.junit.Assert.assertEquals;

import java.util.StringJoiner;

import org.junit.Test;

public class CircosConfigWriterTest
{
    @Test
    public void testCnaAxisPositions() {
        int maxTracks = 42;

        String positionString = CircosConfigWriter.cnaAxisPositions(maxTracks);

        StringJoiner copyNumberString = new StringJoiner(",");
        for (String relativePosition : positionString.split(",")) {
            long copyNumber = Math.round(maxTracks * Double.parseDouble(relativePosition.replace("r", ""))) + 2;
            copyNumberString.add(String.valueOf(copyNumber));
        }

        assertEquals("3,4,5,6,7,8,9,10,20,30,40", copyNumberString.toString());
    }

}
