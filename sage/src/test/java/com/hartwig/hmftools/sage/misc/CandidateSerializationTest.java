package com.hartwig.hmftools.sage.misc;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;

import static org.junit.Assert.assertEquals;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.sage.append.CandidateSerialization;
import com.hartwig.hmftools.sage.candidate_.Candidate_;
import com.hartwig.hmftools.sage.candidate_.IndexedBases_;
import com.hartwig.hmftools.sage.candidate_.ReadContext_;
import com.hartwig.hmftools.sage.common.IndexedBasesTest;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.common.VariantTier;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderVersion;

public class CandidateSerializationTest
{
    private static final VCFCodec CODEC = createTestCodec();

    @NotNull
    private static VCFCodec createTestCodec()
    {
        VCFCodec codec = new VCFCodec();
        VCFHeader header = new VCFHeader(Sets.newHashSet(), Sets.newHashSet("normal", "tumor"));
        codec.setVCFHeader(header, VCFHeaderVersion.VCF4_2);
        return codec;
    }

    @NotNull
    public static Candidate_ decode(String line)
    {
        VariantContext context = CODEC.decode(line);
        IndexedBases_ cheatRefBases = CandidateSerialization.readBases(context);
        return CandidateSerialization.toCandidate(context, cheatRefBases, cheatRefBases);
    }

    @Test
    public void testSerialization()
    {
        final VariantTier expectedTier = VariantTier.HOTSPOT;
        final String expectedRepeat = "AT";
        final int expectedRepeatCount = 2;
        final String expectedMH = "ATGA";
        final int expectedIndex = 1;
        final int expositionPosition = 1000;
        final IndexedBases_ refBases = IndexedBasesTest.createIndexedBases(expositionPosition, expectedIndex, "AA", "TA", "ATG", "CG", "TT");
        final IndexedBases_ readBases = IndexedBasesTest.createIndexedBases(expositionPosition, expectedIndex, "AA", "TA", "ACG", "CG", "TT");

        final ReadContext_ readContext = new ReadContext_(
                expositionPosition, expectedRepeat, expectedRepeatCount, expectedMH, readBases, false);

        final SimpleVariant variant = new SimpleVariant(CHR_1, expositionPosition, "T", "C");

        final Candidate_ candidate = new Candidate_(expectedTier, variant, readContext, 2,0);

        final VariantContext serialized = toContext(candidate);
        final IndexedBases_ deserializedReadBases = CandidateSerialization.readBases(serialized);
        final Candidate_ deserialized = CandidateSerialization.toCandidate(serialized, deserializedReadBases, refBases);

        assertEqual(candidate, deserialized);
        assertEquals(5, candidate.readContext().readBasesPositionIndex());
        assertEquals(3, deserialized.readContext().readBasesPositionIndex());
    }

    private static void assertEqual(Candidate_ expected, Candidate_ victim)
    {
        assertEquals(expected.tier(), victim.tier());
        assertEquals(expected.position(), victim.position());
        assertEquals(expected.chromosome(), victim.chromosome());
        assertEquals(expected.minNumberOfEvents(), victim.minNumberOfEvents());
        assertEquals(expected.readContext().Repeat, victim.readContext().Repeat);
        assertEquals(expected.readContext().RepeatCount, victim.readContext().RepeatCount);
        assertEquals(expected.readContext().microhomology(), victim.readContext().microhomology());
        assertEquals(expected.readContext().leftFlankString(), victim.readContext().leftFlankString());
        assertEquals(expected.readContext().coreString(), victim.readContext().coreString());
        assertEquals(expected.readContext().rightFlankString(), victim.readContext().rightFlankString());
    }

    private static VariantContext toContext(Candidate_ candidate)
    {
        VariantContextBuilder builder = CandidateSerialization.toContext(candidate);

        Genotype genotype = new GenotypeBuilder("SAMPLE").DP(1000).make();
        builder.genotypes(genotype);

        return builder.make();
    }

}
