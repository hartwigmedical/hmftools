package com.hartwig.hmftools.sage.append;

import static com.hartwig.hmftools.common.sage.SageMetaData.TIER;
import static com.hartwig.hmftools.sage.SageConstants.MATCHING_BASE_QUALITY;
import static com.hartwig.hmftools.sage.vcf.VariantVCF.RAW_DEPTH;
import static com.hartwig.hmftools.sage.vcf.VariantVCF.READ_CONTEXT;
import static com.hartwig.hmftools.sage.vcf.VariantVCF.READ_CONTEXT_EVENTS;
import static com.hartwig.hmftools.sage.vcf.VariantVCF.READ_CONTEXT_INDEX;
import static com.hartwig.hmftools.sage.vcf.VariantVCF.READ_CONTEXT_LEFT_FLANK;
import static com.hartwig.hmftools.sage.vcf.VariantVCF.READ_CONTEXT_MICRO_HOMOLOGY;
import static com.hartwig.hmftools.sage.vcf.VariantVCF.READ_CONTEXT_REPEAT_COUNT;
import static com.hartwig.hmftools.sage.vcf.VariantVCF.READ_CONTEXT_REPEAT_SEQUENCE;
import static com.hartwig.hmftools.sage.vcf.VariantVCF.READ_CONTEXT_RIGHT_FLANK;

import java.util.Arrays;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.common.IndexedBases;
import com.hartwig.hmftools.sage.common.ReadContext;
import com.hartwig.hmftools.sage.common.RefSequence;
import com.hartwig.hmftools.sage.common.VariantTier;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public final class CandidateSerialization
{
    public static VariantHotspot toVariantHotspot(final VariantContext context)
    {
        return ImmutableVariantHotspotImpl.builder()
                .chromosome(context.getContig())
                .position(context.getStart())
                .ref(context.getReference().getBaseString())
                .alt(context.getAlternateAllele(0).getBaseString())
                .build();
    }

    public static Candidate toCandidate(final VariantContext context, final RefSequence refGenome)
    {
        final IndexedBases readBases = readBases(context);
        return toCandidate(context, readBases, refGenome.alignment());
    }

    public static IndexedBases readBases(final VariantContext context)
    {
        final int position = context.getStart();

        final String leftFlank = context.getAttributeAsString(READ_CONTEXT_LEFT_FLANK, Strings.EMPTY);
        final String core = context.getAttributeAsString(READ_CONTEXT, Strings.EMPTY);
        final String rightFlank = context.getAttributeAsString(READ_CONTEXT_RIGHT_FLANK, Strings.EMPTY);
        final int readContextIndex = context.getAttributeAsInt(READ_CONTEXT_INDEX, 0);

        final int leftCoreIndex = leftFlank.length();
        final int rightFlankIndex = leftFlank.length() + core.length();
        final int rightCoreIndex = rightFlankIndex - 1;
        final byte[] bases = new byte[leftFlank.length() + core.length() + rightFlank.length()];
        System.arraycopy(leftFlank.getBytes(), 0, bases, 0, leftFlank.length());
        System.arraycopy(core.getBytes(), 0, bases, leftCoreIndex, core.length());
        System.arraycopy(rightFlank.getBytes(), 0, bases, rightFlankIndex, rightFlank.length());

        return new IndexedBases(position,
                readContextIndex + leftFlank.length(),
                leftFlank.length(),
                rightCoreIndex,
                Math.max(leftFlank.length(), rightFlank.length()),
                bases);
    }

    public static Candidate toCandidate(final VariantContext context, final IndexedBases readBases, final IndexedBases refBases)
    {
        final VariantHotspot variant = toVariantHotspot(context);

        final VariantTier tier = VariantTier.valueOf(context.getAttributeAsString(TIER, "LOW_CONFIDENCE"));
        final int repeatCount = context.getAttributeAsInt(READ_CONTEXT_REPEAT_COUNT, 0);
        final String repeat = context.getAttributeAsString(READ_CONTEXT_REPEAT_SEQUENCE, Strings.EMPTY);
        final String mh = context.getAttributeAsString(READ_CONTEXT_MICRO_HOMOLOGY, Strings.EMPTY);

        final ReadContext readContext = new ReadContext(context.getStart(), repeat, repeatCount, mh, readBases, false);

        int maxDepth = 0;
        for(Genotype genotype : context.getGenotypes().immutable())
        {
            maxDepth = Math.max(maxDepth, genotype.getDP());
            maxDepth = Math.max(maxDepth, genotype.getAttributeAsInt(RAW_DEPTH, 0));
        }

        return new Candidate(tier, variant, readContext, maxDepth, context.getAttributeAsInt(READ_CONTEXT_EVENTS, 0), 0);
    }

    public static VariantContextBuilder toContext(final Candidate candidate)
    {
        final List<Allele> alleles = createAlleles(candidate.variant());
        final ReadContext readContext = candidate.readContext();

        final VariantContextBuilder builder = new VariantContextBuilder().chr(candidate.chromosome())
                .source("SAGE")
                .start(candidate.position())
                .attribute(TIER, candidate.tier())
                .attribute(READ_CONTEXT, candidate.readContext().toString())
                .attribute(READ_CONTEXT_LEFT_FLANK, candidate.readContext().leftFlankString())
                .attribute(READ_CONTEXT_RIGHT_FLANK, candidate.readContext().rightFlankString())
                .attribute(READ_CONTEXT_INDEX, readContext.readBasesPositionIndex() - readContext.readBasesLeftCentreIndex())
                .attribute(READ_CONTEXT_EVENTS, candidate.minNumberOfEvents())
                .computeEndFromAlleles(alleles, (int) candidate.position())
                .alleles(alleles);

        if(!readContext.microhomology().isEmpty())
        {
            builder.attribute(READ_CONTEXT_MICRO_HOMOLOGY, readContext.microhomology());
        }

        if(readContext.RepeatCount > 0)
        {
            builder.attribute(READ_CONTEXT_REPEAT_COUNT, readContext.RepeatCount)
                    .attribute(READ_CONTEXT_REPEAT_SEQUENCE, readContext.Repeat);
        }

        return builder;
    }

    @NotNull
    private static List<Allele> createAlleles(final VariantHotspot variant)
    {
        final Allele ref = Allele.create(variant.ref(), true);
        final Allele alt = Allele.create(variant.alt(), false);
        return Lists.newArrayList(ref, alt);
    }
}
