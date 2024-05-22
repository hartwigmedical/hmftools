package com.hartwig.hmftools.sage.evidence;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.sage.common.TestUtils.buildCigarString;
import static com.hartwig.hmftools.sage.common.TestUtils.buildSamRecord;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.sage.common.RefSequence;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.common.VariantReadContextBuilder;
import com.hartwig.hmftools.sage.vcf.CandidateSerialisation;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.sage.common.VariantTier;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class CandidateSerialisationTest
{
    private static MockRefGenome REF_GENOME = new MockRefGenome();

    @Test
    public void testCandidateCreationSnv()
    {
        final String refBases =
            //             10        20        30        40        50
            //   0123456789012345678901234567890123456789012345678901
                "CGCAATATTCGGGTGGGAGTGACCCGATTTACCCGGTGCGTTCGTCACCGCTGTCT";

        RefSequence refSequence = new RefSequence(0, refBases.getBytes());

        REF_GENOME.RefGenomeMap.put(CHR_1, refBases);

        // test 1: a repeat starting at the variant
        String readBases = refBases.substring(0, 21) + "C" + refBases.substring(22, 51);

        String readCigar = buildCigarString(readBases.length());
        SAMRecord read = buildSamRecord(1, readCigar, readBases);

        SimpleVariant var = new SimpleVariant(CHR_1, 21, "A", "C");

        VariantReadContextBuilder builder = new VariantReadContextBuilder(5);
        VariantReadContext readContext = builder.createContext(var, read, 21, refSequence);

        assertTrue(readContext.isValid());

        Candidate candidate = new Candidate(VariantTier.PANEL, readContext, 2, 0);

        VariantContextBuilder variantContextBuilder = CandidateSerialisation.toContext(candidate);

        VariantContext variantContext = variantContextBuilder.make();

        Candidate recreatedCandidate = CandidateSerialisation.toCandidate(variantContext, REF_GENOME);

        VariantReadContext createdContext = recreatedCandidate.readContext();

        assertEquals(createdContext.readCigar(), readContext.readCigar());
        assertEquals(createdContext.AlignmentStart, readContext.AlignmentStart);
        assertEquals(createdContext.AlignmentEnd, readContext.AlignmentEnd);
        assertEquals(createdContext.refBases(), readContext.refBases());
        assertEquals(createdContext.readBases(), readContext.readBases());

        // TODO: fix core position and other things
        // assertEquals(createdContext.CorePositionStart, readContext.CorePositionStart);
        // assertEquals(createdContext.CorePositionEnd, readContext.CorePositionEnd);
    }

    // TODO: add tests for INDELs


    // TODO: add tests for backwards compatibility - ie when no read CIGAR and alignment is available
}
