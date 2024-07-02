package com.hartwig.hmftools.sage.vcf;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

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
import static com.hartwig.hmftools.sage.common.VariantTier.LOW_CONFIDENCE;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_CORE;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_EVENTS;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_INDEX;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_INFO;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_LEFT_FLANK;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_RIGHT_FLANK;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_UPDATED;

import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.common.RefSequence;
import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.common.VariantReadContextBuilder;
import com.hartwig.hmftools.sage.common.VariantTier;

import org.apache.logging.log4j.util.Strings;

import htsjdk.samtools.SAMRecord;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public final class CandidateSerialisation
{
    public static final int PRE_v3_5_FLANK_EXTENSION_LENGTH = 50;

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

        ReadContextVcfInfo readContextVcfInfo = null;
        boolean buildFromOldTags = !context.hasAttribute(READ_CONTEXT_INFO);

        if(buildFromOldTags)
        {
            // support for versions 3.0 -> 3.4
            readContextVcfInfo = buildReadContextFromOldVcfTags(variant, context, refSequence);
        }
        else
        {
            readContextVcfInfo = ReadContextVcfInfo.fromVcfTag(context.getAttributeAsString(READ_CONTEXT_INFO, ""));
        }

        VariantReadContextBuilder builder = new VariantReadContextBuilder(DEFAULT_FLANK_LENGTH);

        SAMRecord record = new SAMRecord(null);
        record.setAlignmentStart(readContextVcfInfo.AlignmentStart);
        record.setCigarString(readContextVcfInfo.Cigar);
        record.setReadBases(readContextVcfInfo.readBases().getBytes());
        record.setReadName("CANDIDATE");

        VariantReadContext readContext = builder.createContext(variant, record, readContextVcfInfo.VarIndex, refSequence);

        if(readContext == null)
        {
            SG_LOGGER.error("variant({}) failed to recreate read context", variant);
            return null;
        }

        // TEMP: tracking of changes
        if(buildFromOldTags)
        {
            // old read base length
            String leftFlank = context.getAttributeAsString(READ_CONTEXT_LEFT_FLANK, Strings.EMPTY);
            String core = context.getAttributeAsString(READ_CONTEXT_CORE, Strings.EMPTY);
            String rightFlank = context.getAttributeAsString(READ_CONTEXT_RIGHT_FLANK, Strings.EMPTY);
            int oldReadBaseLength = leftFlank.length() + core.length() + rightFlank.length();

            if(readContext.totalLength() != oldReadBaseLength)
                context.getCommonInfo().putAttribute(READ_CONTEXT_UPDATED, true);
        }

        return new Candidate(tier, readContext, context.getAttributeAsInt(READ_CONTEXT_EVENTS, 0), 0);
    }

    @VisibleForTesting
    public static ReadContextVcfInfo buildReadContextFromOldVcfTags(
            final SimpleVariant variant, final VariantContext context, final RefSequence refSequence)
    {
        String leftFlank = context.getAttributeAsString(READ_CONTEXT_LEFT_FLANK, Strings.EMPTY);
        String core = context.getAttributeAsString(READ_CONTEXT_CORE, Strings.EMPTY);
        String rightFlank = context.getAttributeAsString(READ_CONTEXT_RIGHT_FLANK, Strings.EMPTY);

        int varCoreIndex = context.getAttributeAsInt(READ_CONTEXT_INDEX, 0);

        // variant index is defined as the variant's index in the core
        int varReadIndex = varCoreIndex + leftFlank.length();

        // eg position = 112, left flank = 10 bases, core = 5 bases, variant core index = 2
        // ie left flank = 100-109, index 0-9, core = 110-114, index 10-14, var read index = 12

        int readBaseLength = leftFlank.length() + core.length() + rightFlank.length();

        int alignmentStart = variant.Position - varReadIndex;
        int alignmentEnd = alignmentStart + readBaseLength - 1 - variant.indelLength();

        // build a buffer around the flanks to allow for any need to expand the core for repeats and homology
        int refBaseBuffer = PRE_v3_5_FLANK_EXTENSION_LENGTH;

        int extendedAlignmentStart = max(alignmentStart - refBaseBuffer, 1);
        String extendedLeft = refSequence.positionBases(extendedAlignmentStart, alignmentStart - 1);
        String extendedRight = refSequence.positionBases(alignmentEnd + 1, alignmentEnd + refBaseBuffer);
        int varIndexOffset = alignmentStart - extendedAlignmentStart;

        varReadIndex += varIndexOffset;
        leftFlank = extendedLeft + leftFlank;
        rightFlank = rightFlank + extendedRight;

        int newReadBaseLength = leftFlank.length() + core.length() + rightFlank.length();

        String readCigar = "";

        if(variant.isInsert())
        {
            // read bases = 10 + 2 + insert of 5 + 2 + 10 = 29, var index = 11, bases after insert = 12
            int preInsertBaseLength = varReadIndex + 1;
            int postInsertBaseLength = newReadBaseLength - varReadIndex - variant.altLength() - 1;
            readCigar = format("%dM%dI%dM", preInsertBaseLength, variant.indelLength(), postInsertBaseLength);
        }
        else if(variant.isDelete())
        {
            // read bases = 10 + 2 + del of 5 + 2 + 10 = 24, var index = 11, bases after insert = 12
            int preInsertBaseLength = varReadIndex + 1;
            int postInsertBaseLength = newReadBaseLength - varReadIndex - 1;
            readCigar = format("%dM%dD%dM", preInsertBaseLength, abs(variant.indelLength()), postInsertBaseLength);
        }
        else
        {
            readCigar = format("%dM", newReadBaseLength);
        }

        return new ReadContextVcfInfo(extendedAlignmentStart, varReadIndex, leftFlank, core, rightFlank, readCigar);
    }
}
