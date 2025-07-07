package com.hartwig.hmftools.common.sv;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.sv.SvUtils.formSvType;
import static com.hartwig.hmftools.common.sv.SvUtils.isIndel;
import static com.hartwig.hmftools.common.sv.SvUtils.isShortLocalDelDupIns;
import static com.hartwig.hmftools.common.sv.SvVcfTags.ALLELE_FRACTION;
import static com.hartwig.hmftools.common.sv.SvVcfTags.ASM_LINKS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.CIPOS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.HOMSEQ;
import static com.hartwig.hmftools.common.sv.SvVcfTags.HOTSPOT;
import static com.hartwig.hmftools.common.sv.SvVcfTags.IHOMPOS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.INFERRED;
import static com.hartwig.hmftools.common.sv.SvVcfTags.INSALN;
import static com.hartwig.hmftools.common.sv.SvVcfTags.MATE_ID;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REF_DEPTH;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REF_DEPTH_PAIR;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SEG_ALIGN_LENGTH;
import static com.hartwig.hmftools.common.sv.SvVcfTags.TOTAL_FRAGS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REPEAT_MASK_REPEAT_CLASS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REPEAT_MASK_COVERAGE;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REPEAT_MASK_ORIENTATION;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REPEAT_MASK_REPEAT_TYPE;
import static com.hartwig.hmftools.common.sv.VariantAltInsertCoords.BREAKEND_REGEX;
import static com.hartwig.hmftools.common.sv.VariantAltInsertCoords.SINGLE_BREAKEND_BYTE;
import static com.hartwig.hmftools.common.sv.VariantAltInsertCoords.SINGLE_BREAKEND_STR;
import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_REV;
import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_FWD;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.getGenotypeAttributeAsDouble;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.getGenotypeAttributeAsInt;

import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.variant.filter.AlwaysPassFilter;
import com.hartwig.hmftools.common.variant.filter.HumanChromosomeFilter;

import org.jetbrains.annotations.Nullable;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.CompoundFilter;
import htsjdk.variant.variantcontext.filter.VariantContextFilter;

public class StructuralVariantFactory
{
    private final Map<String,VariantContext> mUnmatchedVariants;
    private final List<StructuralVariant> mCompleteVariants;

    private final CompoundFilter mFilter;

    private boolean mOrdinalsSet;
    private int mReferenceGenotypeOrdinal;
    private int mTumorGenotypeOrdinal;

    private static final int GENOTYPE_ORDINAL_NONE = -1;

    public static StructuralVariantFactory build(final VariantContextFilter filter)
    {
        CompoundFilter compoundfilter = new CompoundFilter(true);
        compoundfilter.add(new HumanChromosomeFilter());
        compoundfilter.add(filter);
        return new StructuralVariantFactory(compoundfilter);
    }

    public static StructuralVariantFactory build()
    {
        CompoundFilter compoundfilter = new CompoundFilter(true);
        compoundfilter.add(new AlwaysPassFilter());
        return new StructuralVariantFactory(compoundfilter);
    }

    public StructuralVariantFactory(final CompoundFilter filter)
    {
        mFilter = filter;

        mUnmatchedVariants = Maps.newHashMap();
        mCompleteVariants = Lists.newArrayList();
        mReferenceGenotypeOrdinal = GENOTYPE_ORDINAL_NONE;
        mTumorGenotypeOrdinal = GENOTYPE_ORDINAL_NONE;
        mOrdinalsSet = false;
    }

    public void setGenotypeOrdinals(int referenceOrdinal, int tumorOrdinal)
    {
        mOrdinalsSet = true;
        mReferenceGenotypeOrdinal = referenceOrdinal;
        mTumorGenotypeOrdinal = tumorOrdinal;
    }

    public static String mateId(final VariantContext context)
    {
        return context.getAttributeAsString(MATE_ID, null);
    }

    public void clear()
    {
        mCompleteVariants.clear();
        mUnmatchedVariants.clear();
    }

    public static boolean isSingleBreakend(final VariantContext context)
    {
        final byte[] altBases = context.getAlternateAllele(0).getDisplayBases();
        return isSingleBreakend(altBases);
    }

    public static boolean isSingleBreakend(final byte[] altBases)
    {
        return altBases.length > 0 && (altBases[0] == SINGLE_BREAKEND_BYTE || altBases[altBases.length - 1] == SINGLE_BREAKEND_BYTE);
    }

    public static byte parseSingleOrientation(final VariantContext context)
    {
        final String alt = context.getAlternateAllele(0).getDisplayString();
        return alt.startsWith(SINGLE_BREAKEND_STR) ? ORIENT_REV : ORIENT_FWD;
    }

    public static byte parseSvOrientation(final VariantContext context)
    {
        final String alt = context.getAlternateAllele(0).getDisplayString();
        final Matcher match = BREAKEND_REGEX.matcher(alt);

        if(!match.matches())
            return (byte)0;

        return match.group(1).length() > 0 ? ORIENT_FWD : ORIENT_REV;
    }

    public void addVariantContext(final VariantContext context)
    {
        if(mFilter.test(context))
        {
            if(isSingleBreakend(context))
            {
                mCompleteVariants.add(createSingleBreakend(context));
            }
            else
            {
                String mate = context.getAttributeAsString(MATE_ID, null);
                if(mUnmatchedVariants.containsKey(mate))
                {
                    mCompleteVariants.add(createSV(mUnmatchedVariants.remove(mate), context));
                }
                else
                {
                    mUnmatchedVariants.put(context.getID(), context);
                }
            }
        }
    }

    public List<StructuralVariant> results() { return mCompleteVariants; }

    public List<VariantContext> unmatched() { return Lists.newArrayList(mUnmatchedVariants.values()); }

    public StructuralVariant createSV(final VariantContext contextStart, final VariantContext contextEnd)
    {
        final String alt = contextStart.getAlternateAllele(0).getDisplayString();
        final Matcher match = BREAKEND_REGEX.matcher(alt);
        if(!match.matches())
        {
            throw new IllegalArgumentException(String.format("ALT %s is not in breakend notation", alt));
        }

        // Local orientation determined by the positionin of the anchoring bases
        byte startOrientation = (match.group(1).length() > 0 ? ORIENT_FWD : ORIENT_REV);

        // Other orientation determined by the direction of the brackets
        byte endOrientation = (match.group(2).equals("]") ? ORIENT_FWD : ORIENT_REV);

        // Grab the inserted sequence by removing 1 base from the reference anchoring bases
        String insertedSequence = match.group(1).length() > 0 ?
                match.group(1).substring(1) : match.group(4).substring(0, match.group(4).length() - 1);

        StructuralVariantType svType = StructuralVariantType.fromContext(contextStart);

        VariantContext legContextStart = contextStart;
        VariantContext legContextEnd = contextEnd;

        // check for same-base DUP and need to switch context info
        if(contextStart.getContig().equals(contextEnd.getContig()) && startOrientation != endOrientation
                && contextStart.getStart() == contextEnd.getStart() && startOrientation == ORIENT_FWD)
        {
            legContextStart = contextEnd;
            legContextEnd = contextStart;
            startOrientation = ORIENT_REV;
            endOrientation = ORIENT_FWD;
        }

        if(svType == null)
        {
            // infer from attributes
            svType = formSvType(
                    legContextStart.getContig(), legContextEnd.getContig(), legContextStart.getStart(), legContextEnd.getStart(),
                    Orientation.fromByte(startOrientation), Orientation.fromByte(endOrientation), !insertedSequence.isEmpty());
        }

        int indelLength = isIndel(svType) ? abs(contextStart.getStart() - contextEnd.getStart()) : 0;
        boolean isSmallDelDup = isShortLocalDelDupIns(svType, indelLength);

        StructuralVariantLeg startLeg = setLegCommon(legContextStart, isSmallDelDup, startOrientation)
                .position(legContextStart.getStart())
                .homology(legContextStart.getAttributeAsString(HOMSEQ, ""))
                .build();

        StructuralVariantLeg endLeg = setLegCommon(legContextEnd, isSmallDelDup, endOrientation)
                .position(legContextEnd.getStart())
                .homology(legContextEnd.getAttributeAsString(HOMSEQ, ""))
                .build();

        ImmutableStructuralVariantImpl.Builder svBuilder = setCommon(contextStart);

        svBuilder.start(startLeg)
                .end(endLeg)
                .mateId(legContextEnd.getID())
                .insertSequence(insertedSequence)
                .type(svType)
                .filter(filters(legContextStart, legContextEnd))
                .startContext(legContextStart)
                .endContext(legContextEnd);

        svBuilder.startLinkedBy(parseAssemblyLinks(legContextStart));
        svBuilder.endLinkedBy(parseAssemblyLinks(legContextEnd));

        return svBuilder.build();
    }

    private static String parseAssemblyLinks(final VariantContext variantContext)
    {
        return trimStringListValue(variantContext.getAttributeAsString(ASM_LINKS, ""));
    }

    public static String trimStringListValue(final String listValue)
    {
        if(listValue.isEmpty())
            return listValue;

        return listValue.replaceAll("\\[", "").replaceAll("]", "").replaceAll(" ", "");
    }

    public StructuralVariant createSingleBreakend(final VariantContext context)
    {
        final String alt = context.getAlternateAllele(0).getDisplayString();

        // local orientation determined by the positioning of the anchoring bases
        final byte orientation = alt.startsWith(".") ? ORIENT_REV : ORIENT_FWD;
        final int refLength = context.getReference().length();

        final String insertedSequence = orientation == -1 ?
                alt.substring(1, alt.length() - refLength) : alt.substring(refLength, alt.length() - 1);

        final StructuralVariantLeg startLeg = setLegCommon(context, false, orientation)
                .homology("")
                .build();

        return setCommon(context)
                .start(startLeg)
                .insertSequence(insertedSequence)
                .type(context.hasAttribute(INFERRED) ? StructuralVariantType.INF : StructuralVariantType.SGL)
                .filter(filters(context, null))
                .startContext(context)
                .build();
    }

    private ImmutableStructuralVariantImpl.Builder setCommon(final VariantContext context)
    {
        ImmutableStructuralVariantImpl.Builder builder = ImmutableStructuralVariantImpl.builder();

        double qualityScore = context.getPhredScaledQual();

        String insSequenceAlignments = trimStringListValue(context.getAttributeAsString(INSALN, ""));

        builder.id(context.getID())
                .hotspot(context.getAttributeAsBoolean(HOTSPOT, false))
                .event("")
                .qualityScore(qualityScore)
                .insertSequenceAlignments(insSequenceAlignments);

       if(context.hasAttribute(REPEAT_MASK_REPEAT_CLASS))
        {
            builder.insertSequenceRepeatClass(context.getAttributeAsString(REPEAT_MASK_REPEAT_CLASS, ""))
                    .insertSequenceRepeatType(context.getAttributeAsString(REPEAT_MASK_REPEAT_TYPE, ""))
                    .insertSequenceRepeatOrientation(context.getAttributeAsString(REPEAT_MASK_ORIENTATION, "+").equals("+")
                            ? (byte) 1
                            : (byte) -1)
                    .insertSequenceRepeatCoverage(context.getAttributeAsDouble(REPEAT_MASK_COVERAGE, 0));
        }
        return builder;
    }

    private ImmutableStructuralVariantLegImpl.Builder setLegCommon(final VariantContext context, boolean ignoreRefpair, byte orientation)
    {
        ImmutableStructuralVariantLegImpl.Builder builder = ImmutableStructuralVariantLegImpl.builder();
        builder.chromosome(context.getContig());
        builder.position(context.getStart());
        builder.orientation(orientation);

        int ciLeft = 0;
        int ciRight = 0;
        if(context.hasAttribute(CIPOS))
        {
            final List<Integer> cipos = context.getAttributeAsIntList(CIPOS, 0);
            if(cipos.size() == 2)
            {
                ciLeft = cipos.get(0);
                ciRight = cipos.get(1);
            }
        }

        builder.startOffset(ciLeft);
        builder.endOffset(ciRight);

        List<Integer> alignedSegmentLengths = context.getAttributeAsIntList(SEG_ALIGN_LENGTH, 0);
        int maxAnchorLength = alignedSegmentLengths.stream().mapToInt(x -> x.intValue()).max().orElse(0);
        builder.anchoringSupportDistance(maxAnchorLength);

        builder.inexactHomologyOffsetStart(0);
        builder.inexactHomologyOffsetEnd(0);
        if(context.hasAttribute(IHOMPOS))
        {
            final List<Integer> ihompos = context.getAttributeAsIntList(IHOMPOS, 0);
            if(ihompos.size() == 2)
            {
                builder.inexactHomologyOffsetStart(ihompos.get(0));
                builder.inexactHomologyOffsetEnd(ihompos.get(1));
            }
        }

        int referenceOrdinal;
        int tumorOrdinal;

        if(mOrdinalsSet)
        {
            referenceOrdinal = mReferenceGenotypeOrdinal;
            tumorOrdinal = mTumorGenotypeOrdinal;
        }
        else
        {
            referenceOrdinal = context.getGenotypes().size() == 1 ? GENOTYPE_ORDINAL_NONE : 0;
            tumorOrdinal = context.getGenotypes().size() == 1 ? 0 : 1;
        }

        if(referenceOrdinal >= 0 && context.getGenotype(referenceOrdinal) != null)
        {
            Genotype genotype = context.getGenotype(referenceOrdinal);

            int totalFrags = getGenotypeAttributeAsInt(genotype, TOTAL_FRAGS, 0);
            int refFrags = getGenotypeAttributeAsInt(genotype, REF_DEPTH, 0);
            int refPairFrags = getGenotypeAttributeAsInt(genotype, REF_DEPTH_PAIR, 0);
            builder.normalVariantFragmentCount(totalFrags);
            builder.normalReferenceFragmentCount(refFrags + (ignoreRefpair ? 0 : refPairFrags));
        }

        if(tumorOrdinal >= 0 && context.getGenotype(tumorOrdinal) != null)
        {
            Genotype genotype = context.getGenotype(tumorOrdinal);

            int totalFrags = getGenotypeAttributeAsInt(genotype, TOTAL_FRAGS, 0);
            int refFrags = getGenotypeAttributeAsInt(genotype, REF_DEPTH, 0);
            int refPairFrags = getGenotypeAttributeAsInt(genotype, REF_DEPTH_PAIR, 0);
            builder.tumorVariantFragmentCount(totalFrags);
            builder.tumorReferenceFragmentCount(refFrags + (ignoreRefpair ? 0 : refPairFrags));

            double af = getGenotypeAttributeAsDouble(genotype, ALLELE_FRACTION, 0);
            builder.alleleFrequency(af);
        }

        return builder;
    }

    private static String filters(final VariantContext context, @Nullable VariantContext pairedContext)
    {
        final Set<String> filters = new HashSet<>(context.getFilters());
        if(pairedContext != null)
        {
            filters.addAll(pairedContext.getFilters());
        }
        if(filters.size() > 1)
        {
            // Doesn't pass if a filter is applied to either of the two records
            filters.remove(PASS);
        }
        if(filters.size() == 0)
        {
            filters.add(PASS);
        }
        return filters.stream().sorted().collect(Collectors.joining(";"));
    }
}
