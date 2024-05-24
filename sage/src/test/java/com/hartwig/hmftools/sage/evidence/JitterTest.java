package com.hartwig.hmftools.sage.evidence;

import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_FLANK_LENGTH;
import static com.hartwig.hmftools.sage.common.TestUtils.buildSamRecord;
import static com.hartwig.hmftools.sage.common.VariantUtils.createSimpleVariant;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.sage.common.RefSequence;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.sage.common.VariantReadContextBuilder;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class JitterTest
{
    @Test
    public void testReadJitter()
    {
        SimpleVariant variant = createSimpleVariant(30, "A", "T");

        //                           10        20        30        40        50        60        70        80        90
        //                 0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
        String refBases = "ATGTCATTTTAGCGCGCATTCCTTCCTTCCAAAAAAAAAATGCTGGCACACACACATTAGTCGTAGATTAGTCGTAG";
        RefSequence refSequence = new RefSequence(0, refBases.getBytes());

        String readBases = refBases.substring(1, 30) + variant.alt() + refBases.substring(31, 71);
        String readCigar = "70M";
        SAMRecord read = buildSamRecord(1, readCigar, readBases);

        VariantReadContextBuilder builder = new VariantReadContextBuilder(DEFAULT_FLANK_LENGTH);

        VariantReadContext readContext = builder.createContext(variant, read, 29, refSequence);
        assertTrue(readContext.isValid());
        assertEquals(23, readContext.VarIndex);
        assertEquals(3, readContext.AllRepeats.size());

        SAMRecord read1 = buildSamRecord(1, readCigar, readBases);

        JitterMatch jitterMatch = JitterData.checkJitter(readContext, read1, 29);
        assertEquals(JitterMatch.NONE, jitterMatch);

        readBases = refBases.substring(1, 30) + "TTCC" + variant.alt() + refBases.substring(31, 71);
        read1 = buildSamRecord(1, readCigar, readBases);

        jitterMatch = JitterData.checkJitter(readContext, read1, 33);
        assertEquals(JitterMatch.LENGTHENED, jitterMatch);

        readBases = refBases.substring(1, 26) + variant.alt() + refBases.substring(31, 71);
        read1 = buildSamRecord(1, readCigar, readBases);

        jitterMatch = JitterData.checkJitter(readContext, read1, 25);
        assertEquals(JitterMatch.SHORTENED, jitterMatch);

        readBases = refBases.substring(1, 30) + variant.alt() + "A" + refBases.substring(31, 71);
        read1 = buildSamRecord(1, readCigar, readBases);

        jitterMatch = JitterData.checkJitter(readContext, read1, 30);
        assertEquals(JitterMatch.LENGTHENED, jitterMatch);

        readBases = refBases.substring(1, 30) + variant.alt() + refBases.substring(32, 71);
        read1 = buildSamRecord(1, readCigar, readBases);

        jitterMatch = JitterData.checkJitter(readContext, read1, 28);
        assertEquals(JitterMatch.SHORTENED, jitterMatch);
    }
}
