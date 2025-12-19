package com.hartwig.hmftools.sage.vcf;

import static com.hartwig.hmftools.common.redux.BaseQualAdjustment.BASE_QUAL_MINIMUM;
import static com.hartwig.hmftools.common.variant.SageVcfTags.MICROHOMOLOGY;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_MICROHOMOLOGY;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_REPEAT_COUNT;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_REPEAT_SEQUENCE;
import static com.hartwig.hmftools.common.variant.SageVcfTags.REPEAT_COUNT;
import static com.hartwig.hmftools.common.variant.SageVcfTags.REPEAT_SEQUENCE;
import static com.hartwig.hmftools.common.variant.SageVcfTags.TIER;
import static com.hartwig.hmftools.common.variant.SageVcfTags.TRINUCLEOTIDE_CONTEXT;
import static com.hartwig.hmftools.sage.SageCommon.APP_NAME;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_FLANK_LENGTH;
import static com.hartwig.hmftools.common.variant.VariantTier.LOW_CONFIDENCE;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_EVENTS;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_INFO;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.redux.BaseQualAdjustment;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.common.RefSequence;
import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.common.variant.SimpleVariant;
import com.hartwig.hmftools.sage.common.VariantReadContextBuilder;
import com.hartwig.hmftools.common.variant.VariantTier;

import htsjdk.samtools.SAMRecord;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public final class CandidateSerialisation
{
    public static VariantContextBuilder toContext(final Candidate candidate)
    {
        final List<Allele> alleles = createAlleles(candidate.variant());

        final VariantReadContext readContext = candidate.readContext();

        ReadContextVcfInfo readContextVcfInfo = new ReadContextVcfInfo(readContext);

        final VariantContextBuilder builder = new VariantContextBuilder().chr(candidate.chromosome())
                .source(APP_NAME.toUpperCase())
                .start(candidate.position())
                .attribute(TIER, candidate.tier())
                .attribute(READ_CONTEXT_INFO, readContextVcfInfo.toVcfTag())
                .attribute(READ_CONTEXT_EVENTS, candidate.minNumberOfEvents())
                .attribute(TRINUCLEOTIDE_CONTEXT, readContext.trinucleotideStr())
                .computeEndFromAlleles(alleles, candidate.position())
                .alleles(alleles);

        if(readContext.hasHomology())
        {
            builder.attribute(READ_CONTEXT_MICROHOMOLOGY, readContext.homologyBases());
            builder.attribute(MICROHOMOLOGY, readContext.homologyBases()); // not set independently
        }

        if(readContext.MaxRepeat != null)
        {
            builder.attribute(READ_CONTEXT_REPEAT_COUNT, readContext.MaxRepeat.Count);
            builder.attribute(READ_CONTEXT_REPEAT_SEQUENCE, readContext.MaxRepeat.Bases);
        }

        if(readContext.refMaxRepeat() != null)
        {
            builder.attribute(REPEAT_COUNT, readContext.refMaxRepeat().Count);
            builder.attribute(REPEAT_SEQUENCE, readContext.refMaxRepeat().Bases);
        }

        return builder;
    }

    private static List<Allele> createAlleles(final SimpleVariant variant)
    {
        final Allele ref = Allele.create(variant.ref(), true);
        final Allele alt = Allele.create(variant.alt(), false);
        return Lists.newArrayList(ref, alt);
    }

    public static Candidate toCandidate(final VariantContext context, final RefSequence refSequence)
    {
        SimpleVariant variant = new SimpleVariant(
                context.getContig(), context.getStart(),
                context.getReference().getBaseString(), context.getAlternateAllele(0).getBaseString());

        VariantTier tier = VariantTier.valueOf(context.getAttributeAsString(TIER, LOW_CONFIDENCE.toString()));

        ReadContextVcfInfo readContextVcfInfo = ReadContextVcfInfo.fromVcfTag(context.getAttributeAsString(READ_CONTEXT_INFO, ""));

        VariantReadContextBuilder builder = new VariantReadContextBuilder(DEFAULT_FLANK_LENGTH);

        SAMRecord record = new SAMRecord(null);
        record.setAlignmentStart(readContextVcfInfo.AlignmentStart);
        record.setCigarString(readContextVcfInfo.Cigar);
        record.setReadBases(readContextVcfInfo.readBases().getBytes());
        record.setReadName("CANDIDATE");

        byte[] minBaseQuals = new byte[record.getReadBases().length];

        for(int i = 0; i < minBaseQuals.length; ++i)
        {
            minBaseQuals[i] = (byte)40; // maximum across technologies
        }

        record.setBaseQualities(minBaseQuals);

        VariantReadContext readContext = builder.createContext(variant, record, readContextVcfInfo.VarIndex, refSequence);

        if(readContext == null || !readContext.isValid())
        {
            SG_LOGGER.warn("variant({}) failed to recreate read context from VCF info({})", variant, readContextVcfInfo);

            if(readContext == null)
            {
                readContext = new VariantReadContext(
                        variant, readContextVcfInfo.AlignmentStart, 0, null, readContextVcfInfo.readBases().getBytes(),
                        Collections.emptyList(), 0, readContextVcfInfo.VarIndex, 0, null, null,
                        Collections.emptyList(), 0, 0);
            }

            readContext.markInvalid();
        }

        return new Candidate(
                tier, readContext, context.getAttributeAsInt(READ_CONTEXT_EVENTS, 0), 0, 0);
    }
}
