package com.hartwig.hmftools.finding.datamodel;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.junit.Test;

public class DriverInterpretationTest
{
    @Test
    public void sortHighToLow()
    {
        assertEquals(List.of(DriverInterpretation.HIGH, DriverInterpretation.MEDIUM, DriverInterpretation.LOW, DriverInterpretation.UNKNOWN),
                Stream.of(DriverInterpretation.LOW, DriverInterpretation.HIGH, DriverInterpretation.UNKNOWN, DriverInterpretation.MEDIUM)
                        .sorted(DriverInterpretation.highToLow()).collect(Collectors.toList()));
    }
}
