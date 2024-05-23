package com.hartwig.hmftools.sage.vcf;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.variant.SageVcfTags.MICROHOMOLOGY;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_MICROHOMOLOGY;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_REPEAT_COUNT;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_REPEAT_SEQUENCE;
import static com.hartwig.hmftools.common.variant.SageVcfTags.REPEAT_COUNT;
import static com.hartwig.hmftools.common.variant.SageVcfTags.REPEAT_SEQUENCE;
import static com.hartwig.hmftools.common.variant.SageVcfTags.TIER;
import static com.hartwig.hmftools.common.variant.SageVcfTags.TRINUCLEOTIDE_CONTEXT;
import static com.hartwig.hmftools.sage.SageCommon.APP_NAME;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_ALIGNMENT;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_CIGAR;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_CORE;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_EVENTS;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_INDEX;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_LEFT_FLANK;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_RIGHT_FLANK;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.CigarUtils;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.common.Microhomology;
import com.hartwig.hmftools.sage.common.RepeatInfo;
import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.common.VariantTier;

import org.apache.logging.log4j.util.Strings;

import htsjdk.samtools.CigarElement;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public final class CandidateSerialisation
{
    public static VariantContextBuilder toContext(final Candidate candidate)
    {
        final List<Allele> alleles = createAlleles(candidate.variant());

        final VariantReadContext readContext = candidate.readContext();

        final VariantContextBuilder builder = new VariantContextBuilder().chr(candidate.chromosome())
                .source(APP_NAME.toUpperCase())
                .start(candidate.position())
                .attribute(TIER, candidate.tier())
                .attribute(READ_CONTEXT_CORE, readContext.coreStr())
                .attribute(READ_CONTEXT_LEFT_FLANK, readContext.leftFlankStr())
                .attribute(READ_CONTEXT_RIGHT_FLANK, readContext.rightFlankStr())
                .attribute(READ_CONTEXT_INDEX, readContext.VarReadIndex)
                .attribute(READ_CONTEXT_CIGAR, readContext.readCigar())
                .attribute(READ_CONTEXT_ALIGNMENT, new int[] { readContext.AlignmentStart, readContext.AlignmentEnd } )
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

    public static Candidate toCandidate(final VariantContext context, final RefGenomeInterface refGenome)
    {
        SimpleVariant variant = new SimpleVariant(
                context.getContig(), context.getStart(),
                context.getReference().getBaseString(), context.getAlternateAllele(0).getBaseString());

        VariantTier tier = VariantTier.valueOf(context.getAttributeAsString(TIER, "LOW_CONFIDENCE"));

        int repeatCount = context.getAttributeAsInt(READ_CONTEXT_REPEAT_COUNT, 0);
        String repeatSequence = context.getAttributeAsString(READ_CONTEXT_REPEAT_SEQUENCE, Strings.EMPTY);
        String homologyStr = context.getAttributeAsString(READ_CONTEXT_MICROHOMOLOGY, Strings.EMPTY);

        String leftFlank = context.getAttributeAsString(READ_CONTEXT_LEFT_FLANK, Strings.EMPTY);
        String core = context.getAttributeAsString(READ_CONTEXT_CORE, Strings.EMPTY);
        String rightFlank = context.getAttributeAsString(READ_CONTEXT_RIGHT_FLANK, Strings.EMPTY);
        int varReadIndex = context.getAttributeAsInt(READ_CONTEXT_INDEX, 0);

        String readCigarStr = context.getAttributeAsString(READ_CONTEXT_CIGAR, Strings.EMPTY);
        List<Integer> readAlignment = context.getAttributeAsIntList(READ_CONTEXT_ALIGNMENT, 0);

        int readAlignmentStart = 0;
        int readAlignmentEnd = 0;

        final String readBases = leftFlank + core + rightFlank;

        List<CigarElement> readCigar;

        if(!readCigarStr.isEmpty() && readAlignment.size() == 2)
        {
            readCigar = CigarUtils.cigarElementsFromStr(readCigarStr);
            readAlignmentStart = readAlignment.get(0);
            readAlignmentEnd = readAlignment.get(1);
        }
        else
        {
            // pre v3.5 the var index was only with reference to the CORE, not including the flanks
            varReadIndex += leftFlank.length();
            readCigar = recreateReadCigar(variant, varReadIndex, readBases.length());
            readAlignmentStart = variant.Position - varReadIndex;
            readAlignmentEnd = readAlignmentStart + readBases.length() - 1 - variant.indelLength();
        }

        int coreIndexStart = leftFlank.length();
        int coreIndexEnd = coreIndexStart + core.length() - 1;

        int leftCoreLength = varReadIndex - coreIndexStart;
        int rightCoreLength = coreIndexEnd - varReadIndex;

        // now should just use core positions
        int[] refPosCoords = new int[] {0, 0}; // getRefBaseCoordinates(variant, core.length(), leftCoreLength, rightCoreLength);

        final byte[] refBases = refGenome.getBases(variant.Chromosome, refPosCoords[SE_START], refPosCoords[SE_END]);

        Microhomology homology = !homologyStr.isEmpty() ? new Microhomology(homologyStr, homologyStr.length()) : null;
        RepeatInfo maxRepeat = !repeatSequence.isEmpty() ? new RepeatInfo(0, repeatSequence, repeatCount) : null;

        // TODO: all repeats, core positions etc etc

        VariantReadContext readContext = new VariantReadContext(
                variant, readAlignmentStart, readAlignmentEnd, refBases, readBases.getBytes(), readCigar,
                coreIndexStart, varReadIndex, coreIndexEnd, homology, maxRepeat, Lists.newArrayList(),
                0, 0, 0, 0);

        return new Candidate(tier, readContext, context.getAttributeAsInt(READ_CONTEXT_EVENTS, 0), 0);
    }

    public static List<CigarElement> recreateReadCigar(final SimpleVariant variant, int varReadIndex, int totalReadLength)
    {
        List<CigarElement> readCigar = Lists.newArrayList();

        if(variant.isIndel())
        {
            readCigar.add(new CigarElement(varReadIndex - 1, M));

            if(variant.isInsert())
                readCigar.add(new CigarElement(variant.Alt.length() - 1, I));
            else
                readCigar.add(new CigarElement(variant.Ref.length() - 1, D));

            int remainingBases = totalReadLength - varReadIndex - variant.Alt.length();
            readCigar.add(new CigarElement(remainingBases, M));
        }
        else
        {
            readCigar.add(new CigarElement(totalReadLength, M));
        }

        return readCigar;
    }
}
