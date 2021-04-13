package com.hartwig.hmftools.sage.candidate;

import static org.junit.Assert.assertEquals;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.read.IndexedBases;
import com.hartwig.hmftools.sage.read.IndexedBasesTest;
import com.hartwig.hmftools.sage.read.ReadContext;
import com.hartwig.hmftools.sage.variant.SageVariantTier;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderVersion;

public class CandidateSerializationTest {

    private static final VCFCodec CODEC = createTestCodec();

    @NotNull
    private static VCFCodec createTestCodec() {
        VCFCodec codec = new VCFCodec();
        VCFHeader header = new VCFHeader(Sets.newHashSet(), Sets.newHashSet("normal", "tumor"));
        codec.setVCFHeader(header, VCFHeaderVersion.VCF4_2);
        return codec;
    }

    @NotNull
    public static Candidate decode(String line) {
        VariantContext context = CODEC.decode(line);
        IndexedBases cheatRefBases = CandidateSerialization.readBases(context);
        return CandidateSerialization.toCandidate(context, cheatRefBases, cheatRefBases);
    }

    @Test
    public void testSerialization() {
        final SageVariantTier expectedTier = SageVariantTier.HOTSPOT;
        final String expectedRepeat = "AT";
        final int expectedRepeatCount = 2;
        final String expectedMH = "ATGA";
        final int expectedIndex = 1;
        final int expositionPosition = 1000;
        final VariantHotspot variant =
                ImmutableVariantHotspotImpl.builder().position(expositionPosition).chromosome("1").ref("T").alt("C").build();
        final IndexedBases refBases = IndexedBasesTest.create(expositionPosition, expectedIndex, "AA", "TA", "ATG", "CG", "TT");
        final IndexedBases readBases = IndexedBasesTest.create(expositionPosition, expectedIndex, "AA", "TA", "ACG", "CG", "TT");
        final ReadContext readContext = new ReadContext(refBases, readBases, expectedRepeatCount, expectedRepeat, expectedMH);

        final Candidate candidate = new Candidate(expectedTier, variant, readContext, 1000, 2);

        final VariantContext serialized = toContext(candidate);
        final IndexedBases deserializedReadBases = CandidateSerialization.readBases(serialized);
        final Candidate deserialized = CandidateSerialization.toCandidate(serialized, deserializedReadBases, refBases);

        assertEqual(candidate, deserialized);
        assertEquals(5, candidate.readContext().readBasesPositionIndex());
        assertEquals(3, deserialized.readContext().readBasesPositionIndex());
    }

    private static void assertEqual(Candidate expected, Candidate victim) {
        assertEquals(expected.tier(), victim.tier());
        assertEquals(expected.position(), victim.position());
        assertEquals(expected.chromosome(), victim.chromosome());
        assertEquals(expected.maxReadDepth(), victim.maxReadDepth());
        assertEquals(expected.minNumberOfEvents(), victim.minNumberOfEvents());
        assertEquals(expected.readContext().repeat(), victim.readContext().repeat());
        assertEquals(expected.readContext().repeatCount(), victim.readContext().repeatCount());
        assertEquals(expected.readContext().microhomology(), victim.readContext().microhomology());
        assertEquals(expected.readContext().leftFlankString(), victim.readContext().leftFlankString());
        assertEquals(expected.readContext().centerBases(), victim.readContext().centerBases());
        assertEquals(expected.readContext().rightFlankString(), victim.readContext().rightFlankString());
    }

    private static VariantContext toContext(Candidate candidate) {
        VariantContextBuilder builder = CandidateSerialization.toContext(candidate);

        Genotype genotype = new GenotypeBuilder("SAMPLE").DP(candidate.maxReadDepth()).make();
        builder.genotypes(genotype);

        return builder.make();
    }

}
