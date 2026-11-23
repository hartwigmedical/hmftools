package com.hartwig.hmftools.cup;

import static java.util.function.Function.identity;
import static java.util.stream.Collectors.toMap;

import static org.junit.Assert.assertEquals;

import java.util.Map;

import com.hartwig.hmftools.cup.drivers.ViralInsertionType;

import org.junit.Test;

public class ViralInsertionTypeTest
{
    // Names of reportable viruses in the Viridian reference
    @Test
    public void testFromVirusName()
    {
        Map<String, ViralInsertionType> expected = Map.of(
                "Human papillomavirus type 16", ViralInsertionType.HPV,
                "Hepatitis B virus", ViralInsertionType.HBV,
                "Merkel cell polyomavirus", ViralInsertionType.MERKEL,
                "Human gammaherpesvirus 4", ViralInsertionType.HERPES,
                "Human gammaherpesvirus 8", ViralInsertionType.HERPES);

        Map<String, ViralInsertionType> actual = expected.keySet().stream()
                .collect(toMap(identity(), ViralInsertionType::fromVirusName));

        assertEquals(expected, actual);
    }
}
