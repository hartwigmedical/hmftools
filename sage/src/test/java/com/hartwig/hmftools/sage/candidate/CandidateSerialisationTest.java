package com.hartwig.hmftools.sage.candidate;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_FLANK_LENGTH;
import static com.hartwig.hmftools.sage.common.TestUtils.buildCigarString;
import static com.hartwig.hmftools.sage.common.TestUtils.buildSamRecord;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_CORE;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_INDEX;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_INFO;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_LEFT_FLANK;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_RIGHT_FLANK;

import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.sage.common.RefSequence;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.common.VariantReadContextBuilder;
import com.hartwig.hmftools.sage.vcf.CandidateSerialisation;
import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.sage.common.VariantTier;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class CandidateSerialisationTest
{
    private static MockRefGenome REF_GENOME = new MockRefGenome();

    private static final String REF_BASES =
        //             10        20        30        40        50
        //   0123456789012345678901234567890123456789012345678901
            "CGCAATATTCGGGTGGGAGTGACCCGATTTACCCGGTGCGTTCGTCACCGCTGTCT";

    private static RefSequence REF_SEQUENCE = new RefSequence(1, REF_BASES.getBytes());

    static
    {
        REF_GENOME.RefGenomeMap.put(CHR_1, REF_BASES);
    }

    @Test
    public void testCandidateCreation()
    {
        // test 1: a simple SNV
        int position = 21;
        String readBases = REF_BASES.substring(1, position) + "C" + REF_BASES.substring(position + 1, 50);

        String readCigar = buildCigarString(readBases.length());
        SAMRecord read = buildSamRecord(1, readCigar, readBases);

        SimpleVariant var = new SimpleVariant(CHR_1, position, "A", "C");

        VariantReadContextBuilder builder = new VariantReadContextBuilder(DEFAULT_FLANK_LENGTH);
        VariantReadContext readContext = builder.createContext(var, read, position - 1, REF_SEQUENCE);

        assertTrue(readContext.isValid());

        Candidate candidate = new Candidate(VariantTier.PANEL, readContext, 2, 0);

        VariantContextBuilder variantContextBuilder = CandidateSerialisation.toContext(candidate);

        VariantContext variantContext = variantContextBuilder.make();

        Candidate recreatedCandidate = CandidateSerialisation.toCandidate(variantContext, REF_SEQUENCE);
        assertNotNull(recreatedCandidate);

        VariantReadContext recreatedContext = recreatedCandidate.readContext();

        assertTrue(recreatedContext.matches(readContext));

        testCandidateCreationOldVersions(readContext, candidate);

        // a delete
        position = 17;
        var = new SimpleVariant(CHR_1, position, "ACCC", "A");

        readBases = REF_BASES.substring(1, position + 1) + REF_BASES.substring(position + 4, 50);
        readCigar = "17M3D29M";
        read = buildSamRecord(1, readCigar, readBases);

        readContext = builder.createContext(var, read, position - 1, REF_SEQUENCE);

        assertTrue(readContext.isValid());

        candidate = new Candidate(VariantTier.PANEL, readContext, 2, 0);

        variantContextBuilder = CandidateSerialisation.toContext(candidate);

        variantContext = variantContextBuilder.make();

        recreatedCandidate = CandidateSerialisation.toCandidate(variantContext, REF_SEQUENCE);
        assertNotNull(recreatedCandidate);

        recreatedContext = recreatedCandidate.readContext();
        assertTrue(recreatedContext.matches(readContext));

        testCandidateCreationOldVersions(readContext, candidate);

        // an insert
        position = 21;
        var = new SimpleVariant(CHR_1, position, "A", "ACCC");

        readBases = REF_BASES.substring(1, position + 1); // up to and including the ref base
        readBases += "CCC" + REF_BASES.substring(position + 2, 50);
        readCigar = "21M3I27M";
        read = buildSamRecord(1, readCigar, readBases);

        readContext = builder.createContext(var, read, position - 1, REF_SEQUENCE);

        assertTrue(readContext.isValid());

        candidate = new Candidate(VariantTier.PANEL, readContext, 2, 0);

        variantContextBuilder = CandidateSerialisation.toContext(candidate);

        variantContext = variantContextBuilder.make();

        recreatedCandidate = CandidateSerialisation.toCandidate(variantContext, REF_SEQUENCE);
        assertNotNull(recreatedCandidate);

        recreatedContext = recreatedCandidate.readContext();
        assertTrue(recreatedContext.matches(readContext));

        testCandidateCreationOldVersions(readContext, candidate);
    }

    private void testCandidateCreationOldVersions(final VariantReadContext readContext, final Candidate candidate)
    {
        VariantContextBuilder variantContextBuilder = CandidateSerialisation.toContext(candidate);

        variantContextBuilder.getAttributes().remove(READ_CONTEXT_INFO);

        variantContextBuilder.getAttributes().put(READ_CONTEXT_LEFT_FLANK, readContext.leftFlankStr());
        variantContextBuilder.getAttributes().put(READ_CONTEXT_CORE, readContext.coreStr());
        variantContextBuilder.getAttributes().put(READ_CONTEXT_RIGHT_FLANK, readContext.rightFlankStr());

        int varIndexInCore = readContext.VarIndex - readContext.leftFlankLength();
        variantContextBuilder.getAttributes().put(READ_CONTEXT_INDEX, varIndexInCore);

        VariantContext variantContext = variantContextBuilder.make();

        Candidate recreatedCandidate = CandidateSerialisation.toCandidate(variantContext, REF_SEQUENCE);
        assertNotNull(recreatedCandidate);

        VariantReadContext recreatedContext = recreatedCandidate.readContext();

        assertTrue(recreatedContext.matches(readContext));
    }
}
