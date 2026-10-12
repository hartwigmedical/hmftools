package com.hartwig.hmftools.viridian.app;

import static com.hartwig.hmftools.viridian.app.ViridianConfig.parseExtraVirusDecoyContigs;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThrows;

import java.util.Set;

import com.hartwig.hmftools.viridian.common.UserInputError;

import org.junit.Test;

public class ViridianConfigTest
{
    @Test
    public void testParseExtraVirusDecoyContigsUnsetIsEmpty()
    {
        assertEquals(Set.of(), parseExtraVirusDecoyContigs(null));
    }

    @Test
    public void testParseExtraVirusDecoyContigsSplitsOnSemicolon()
    {
        assertEquals(Set.of("chrHPV16"), parseExtraVirusDecoyContigs("chrHPV16"));
        assertEquals(Set.of("chrHPV16", "chrHBV"), parseExtraVirusDecoyContigs("chrHPV16;chrHBV"));
    }

    @Test
    public void testParseExtraVirusDecoyContigsBlankContigRejected()
    {
        assertThrows(UserInputError.class, () -> parseExtraVirusDecoyContigs(""));
        assertThrows(UserInputError.class, () -> parseExtraVirusDecoyContigs("chrHPV16;;chrHBV"));
        assertThrows(UserInputError.class, () -> parseExtraVirusDecoyContigs("chrHPV16;"));
        assertThrows(UserInputError.class, () -> parseExtraVirusDecoyContigs("chrHPV16; "));
    }
}
