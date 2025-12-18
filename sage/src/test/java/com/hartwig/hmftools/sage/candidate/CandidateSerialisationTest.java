package com.hartwig.hmftools.sage.candidate;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.variant.SageVcfTags.REPEAT_COUNT;
import static com.hartwig.hmftools.common.variant.SageVcfTags.REPEAT_SEQUENCE;
import static com.hartwig.hmftools.common.variant.SageVcfTags.TIER;
import static com.hartwig.hmftools.common.variant.SageVcfTags.TRINUCLEOTIDE_CONTEXT;
import static com.hartwig.hmftools.sage.SageCommon.APP_NAME;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_FLANK_LENGTH;
import static com.hartwig.hmftools.sage.common.TestUtils.TEST_SAMPLE;
import static com.hartwig.hmftools.sage.common.TestUtils.buildCigarString;
import static com.hartwig.hmftools.sage.common.TestUtils.buildSamRecord;
import static com.hartwig.hmftools.sage.vcf.CandidateSerialisation.toCandidate;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_EVENTS;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_INFO;

import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import static htsjdk.variant.vcf.VCFConstants.ALLELE_FREQUENCY_KEY;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.sage.common.RefSequence;
import com.hartwig.hmftools.common.variant.SimpleVariant;
import com.hartwig.hmftools.sage.common.VariantReadContextBuilder;
import com.hartwig.hmftools.sage.vcf.CandidateSerialisation;
import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.sage.vcf.ReadContextVcfInfo;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class CandidateSerialisationTest
{
    private static final String REF_BASES =
        //             10        20        30        40        50
        //   0123456789012345678901234567890123456789012345678901
            "CGCAATATTCGGGTGGGAGTGACCCGATTTACCCGGTGCGTTCGTCACCGCTGTCT";

    private static RefSequence REF_SEQUENCE = new RefSequence(1, REF_BASES.getBytes());

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

        Candidate recreatedCandidate = toCandidate(variantContext, REF_SEQUENCE);
        assertNotNull(recreatedCandidate);

        VariantReadContext recreatedContext = recreatedCandidate.readContext();

        assertTrue(recreatedContext.matches(readContext));

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

        recreatedCandidate = toCandidate(variantContext, REF_SEQUENCE);
        assertNotNull(recreatedCandidate);

        recreatedContext = recreatedCandidate.readContext();
        assertTrue(recreatedContext.matches(readContext));

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

        recreatedCandidate = toCandidate(variantContext, REF_SEQUENCE);
        assertNotNull(recreatedCandidate);

        recreatedContext = recreatedCandidate.readContext();
        assertTrue(recreatedContext.matches(readContext));
    }

    @Test
    public void testCreateFromVariantContext()
    {
        // public static Candidate toCandidate(final VariantContext context, final RefSequence refSequence)
        GenotypeBuilder genotypeBuilder = new GenotypeBuilder(TEST_SAMPLE);
        genotypeBuilder.AD(new int[] { 10, 10 });
        genotypeBuilder.DP(100);
        genotypeBuilder.alleles(Lists.newArrayList(Allele.NO_CALL, Allele.NO_CALL));
        genotypeBuilder.attribute(ALLELE_FREQUENCY_KEY, 0.1);
        Genotype genotype = genotypeBuilder.make();

        SimpleVariant variant = new SimpleVariant(CHR_1, 17, "AAAAAT", "A");

        List<Allele> alleles = List.of(Allele.create(variant.Ref, true), Allele.create(variant.Alt, false));

        VariantContextBuilder builder = new VariantContextBuilder();
        builder.chr(variant.Chromosome).start(variant.Position);
        builder.alleles(alleles);
        builder.computeEndFromAlleles(alleles, variant.Position);

        builder.genotypes(List.of(genotype));

        // index                                   0123456789 012345678901234567890 1234567890
        builder.attribute(READ_CONTEXT_INFO, "1-13-TTCAGGTTTG-GAAAAAAATATATATATATAA-GTGCATGAAA-17M5D24M");

        VariantContext variantContext = builder.make();

        Candidate candidate = toCandidate(variantContext, REF_SEQUENCE);
        assertNotNull(candidate);
        assertNotNull(candidate.readContext().toString());
    }
}
