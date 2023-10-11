package com.hartwig.hmftools.common.sv;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;

import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import com.google.common.base.Preconditions;
import com.google.common.base.Strings;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.filter.ExcludeCNVFilter;
import com.hartwig.hmftools.common.variant.filter.HumanChromosomeFilter;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.CompoundFilter;
import htsjdk.variant.variantcontext.filter.VariantContextFilter;

public class StructuralVariantFactory
{
    public static final String RECOVERED = "RECOVERED";
    public static final String RECOVERY_METHOD = "RECOVERY_METHOD";
    public static final String RECOVERY_FILTER = "RECOVERY_FILTER";
    public static final String PASS = "PASS";
    public static final String INFERRED = "INFERRED";
    public static final String CIPOS = "CIPOS";
    public static final String SVTYPE = "SVTYPE";
    public static final String PON_FILTER_PON = "PON";
    public static final String IMPRECISE = "IMPRECISE";
    public static final String MATE_ID = "MATEID";
    public static final String PAR_ID = "PARID";

    // set by Gripss
    public static final String TAF = "TAF";
    public static final String HOTSPOT = "HOTSPOT";
    public static final String LOCAL_LINKED_BY = "LOCAL_LINKED_BY";
    public static final String REMOTE_LINKED_BY = "REMOTE_LINKED_BY";
    public static final String LINKED_BY_DELIM = ",";
    public static final String PON_COUNT = "PON_COUNT";

    private static final String INS_SEQ = "SVINSSEQ";
    private static final String HOM_SEQ = "HOMSEQ";
    private static final String SOMATIC_SCORE = "SOMATICSCORE"; // only applicable for Manta and will be removed when fully on GRIDSS
    public static final String IHOMPOS = "IHOMPOS";

    public static final String ALLELE_FRACTION = "AF";
    public static final String ALLELE_FRACTION_DESC = "Allele fraction";

    public static final String SV_FRAGMENT_COUNT = "VF";
    public static final String SGL_FRAGMENT_COUNT = "BVF";

    public static final String REF_READ_COVERAGE = "REF";
    public static final String REF_READ_COVERAGE_DESC = "Count of reads mapping across this breakend";

    public static final String REF_READPAIR_COVERAGE = "REFPAIR";

    public static final String REF_READPAIR_COVERAGE_DESC =
            "Count of reference read pairs spanning this breakend supporting the reference allele";

    private static final String EVENT = "EVENT";

    public static final String UNTEMPLATED_SEQUENCE_ALIGNMENTS = "BEALN";
    public static final String UNTEMPLATED_SEQUENCE_REPEAT_CLASS = "INSRMRC";
    public static final String UNTEMPLATED_SEQUENCE_REPEAT_TYPE = "INSRMRT";
    public static final String UNTEMPLATED_SEQUENCE_REPEAT_ORIENTATION = "INSRMRO";
    public static final String UNTEMPLATED_SEQUENCE_REPEAT_COVERAGE = "INSRMP";

    private static final String ANCHOR_SUPPORT_CIGAR = "SC";

    // set by Purple
    public static final String REF_CONTEXT_FLAG = "REFG";

    // Must match the small deldup threshold in Gripss
    private static final int SMALL_DELDUP_SIZE = 1000;

    public static final Pattern BREAKEND_REGEX = Pattern.compile("^(.*)([\\[\\]])(.+)[\\[\\]](.*)$");
    public static final Pattern SINGLE_BREAKEND_REGEX = Pattern.compile("^(([.].*)|(.*[.]))$");
    private static final byte SINGLE_BREAKEND_BYTE = 46;

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
        compoundfilter.add(new ExcludeCNVFilter());
        compoundfilter.add(filter);
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

    public static Direction parseSingleOrientation(final VariantContext context)
    {
        final String alt = context.getAlternateAllele(0).getDisplayString();
        return alt.startsWith(".") ? Direction.REVERSE : Direction.FORWARDS;
    }

    public static Direction parseSvOrientation(final VariantContext context)
    {
        final String alt = context.getAlternateAllele(0).getDisplayString();
        final Matcher match = BREAKEND_REGEX.matcher(alt);

        if(!match.matches())
            return Direction.UNKNOWN;

        return match.group(1).length() > 0 ? Direction.FORWARDS : Direction.REVERSE;
    }

    public void addVariantContext(final VariantContext context)
    {
        if(mFilter.test(context))
        {
            final StructuralVariantType type = type(context);
            if(type.equals(BND))
            {
                if(isSingleBreakend(context))
                {
                    mCompleteVariants.add(createSingleBreakend(context));
                }
                else
                {
                    String mate = mateId(context);
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
            else
            {
                // not from Gridss and not supported
            }
        }
    }

    @Nullable
    public static String mateId(final VariantContext context)
    {
        String mate = context.getAttributeAsString(MATE_ID, null);

        if(mate == null)
            return context.getAttributeAsString(PAR_ID, null);

        return mate;
    }

    public List<StructuralVariant> results() { return mCompleteVariants; }

    public List<VariantContext> unmatched() { return Lists.newArrayList(mUnmatchedVariants.values()); }

    public void removeUnmatchedVariant(final String id) { mUnmatchedVariants.remove(id); }
    public boolean hasUnmatchedVariant(final String id) { return mUnmatchedVariants.containsKey(id); }

    public StructuralVariant createSV(final VariantContext first, final VariantContext second)
    {
        Preconditions.checkArgument(BND.equals(type(first)));
        Preconditions.checkArgument(BND.equals(type(second)));

        final int start = first.getStart();
        final int end = second.getStart();

        final String alt = first.getAlternateAllele(0).getDisplayString();
        final Matcher match = BREAKEND_REGEX.matcher(alt);
        if(!match.matches())
        {
            throw new IllegalArgumentException(String.format("ALT %s is not in breakend notation", alt));
        }

        // Local orientation determined by the positioning of the anchoring bases
        final byte startOrientation = (byte) (match.group(1).length() > 0 ? POS_ORIENT : NEG_ORIENT);

        // Other orientation determined by the direction of the brackets
        final byte endOrientation = (byte) (match.group(2).equals("]") ? POS_ORIENT : NEG_ORIENT);

        // Grab the inserted sequence by removing 1 base from the reference anchoring bases
        String insertedSequence = match.group(1).length() > 0 ?
                match.group(1).substring(1) : match.group(4).substring(0, match.group(4).length() - 1);

        if(Strings.isNullOrEmpty(insertedSequence))
        {
            insertedSequence = first.getAttributeAsString(INS_SEQ, "");
        }

        final boolean isSmallDelDup = first.getContig().equals(second.getContig())
                && Math.abs(first.getStart() - second.getStart()) <= SMALL_DELDUP_SIZE
                && startOrientation != endOrientation;

        final StructuralVariantLeg startLeg = setLegCommon(first, isSmallDelDup, startOrientation)
                .position(start)
                .homology(first.getAttributeAsString(HOM_SEQ, ""))
                .alleleFrequency(unadjustedAllelicFrequency(first))
                .build();

        final StructuralVariantLeg endLeg = setLegCommon(second, isSmallDelDup, endOrientation)
                .position(end)
                .homology(second.getAttributeAsString(HOM_SEQ, ""))
                .alleleFrequency(unadjustedAllelicFrequency(second))
                .build();

        StructuralVariantType inferredType = BND;
        if(startLeg.chromosome().equals(endLeg.chromosome()))
        {
            if(startLeg.orientation() == endLeg.orientation())
            {
                inferredType = StructuralVariantType.INV;
            }
            else if(startLeg.orientation() == -1)
            {
                inferredType = StructuralVariantType.DUP;
            }
            else if(insertedSequence != null && insertedSequence.length() > 0
                    && Math.abs(endLeg.position() - startLeg.position()) <= 1)
            {
                inferredType = StructuralVariantType.INS;
            }
            else
            {
                inferredType = StructuralVariantType.DEL;
            }
        }

        return setCommon(first)
                .start(startLeg)
                .end(endLeg)
                .mateId(second.getID())
                .insertSequence(insertedSequence)
                .type(inferredType)
                .filter(filters(first, second))
                .startContext(first)
                .endContext(second)
                .build();
    }

    public StructuralVariant createSingleBreakend(final VariantContext context)
    {
        final List<Double> af = context.hasAttribute(TAF) ? context.getAttributeAsDoubleList(TAF, 0.0) : Collections.emptyList();

        final String alt = context.getAlternateAllele(0).getDisplayString();

        // local orientation determined by the positioning of the anchoring bases
        final byte orientation = alt.startsWith(".") ? NEG_ORIENT : POS_ORIENT;
        final int refLength = context.getReference().length();

        final String insertedSequence = orientation == -1 ?
                alt.substring(1, alt.length() - refLength) : alt.substring(refLength, alt.length() - 1);

        final StructuralVariantLeg startLeg = setLegCommon(context, false, orientation)
                .homology("")
                .alleleFrequency(af.size() >= 1 ? af.get(0) : null)
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

        // backwards compatibility for Manta until fully over to GRIDSS
        int somaticScore = context.getAttributeAsInt(SOMATIC_SCORE, 0);
        double qualityScore = context.getPhredScaledQual();

        if(somaticScore > 0)
        {
            qualityScore = somaticScore;
        }

        builder.id(context.getID())
                .recovered(context.getAttributeAsBoolean(RECOVERED, false))
                .hotspot(context.getAttributeAsBoolean(HOTSPOT, false))
                .recoveryMethod(context.getAttributeAsString(RECOVERY_METHOD, null))
                .recoveryFilter(context.getAttributeAsStringList(RECOVERY_FILTER, "").stream().collect(Collectors.joining(",")))
                .event(context.getAttributeAsString(EVENT, null))
                .startLinkedBy(context.getAttributeAsStringList(LOCAL_LINKED_BY, "")
                        .stream()
                        .filter(s -> !Strings.isNullOrEmpty(s))
                        .collect(Collectors.joining(",")))
                .endLinkedBy(context.getAttributeAsStringList(REMOTE_LINKED_BY, "")
                        .stream()
                        .filter(s -> !Strings.isNullOrEmpty(s))
                        .collect(Collectors.joining(",")))
                .imprecise(imprecise(context))
                .qualityScore(qualityScore)
                .insertSequenceAlignments(context.getAttributeAsStringList(UNTEMPLATED_SEQUENCE_ALIGNMENTS, "")
                        .stream()
                        .filter(s -> !Strings.isNullOrEmpty(s))
                        .collect(Collectors.joining(",")));
        if(context.hasAttribute(UNTEMPLATED_SEQUENCE_REPEAT_CLASS))
        {
            builder.insertSequenceRepeatClass(context.getAttributeAsString(UNTEMPLATED_SEQUENCE_REPEAT_CLASS, ""))
                    .insertSequenceRepeatType(context.getAttributeAsString(UNTEMPLATED_SEQUENCE_REPEAT_TYPE, ""))
                    .insertSequenceRepeatOrientation(context.getAttributeAsString(UNTEMPLATED_SEQUENCE_REPEAT_ORIENTATION, "+").equals("+")
                            ? (byte) 1
                            : (byte) -1)
                    .insertSequenceRepeatCoverage(context.getAttributeAsDouble(UNTEMPLATED_SEQUENCE_REPEAT_COVERAGE, 0));
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

        int supportWidth = 0;
        if(context.hasAttribute(ANCHOR_SUPPORT_CIGAR))
        {
            Cigar cigar = TextCigarCodec.decode(context.getAttributeAsString(ANCHOR_SUPPORT_CIGAR, ""));
            for(CigarElement ce : cigar.getCigarElements())
            {
                switch(ce.getOperator())
                {
                    case X:
                    case N:
                        break;
                    case D:
                    case M:
                    case EQ:
                        supportWidth += ce.getLength();
                        break;
                    default:
                        throw new RuntimeException(String.format("Unsupported support interval cigar operator \"%s\"",
                                ce.getOperator().toString()));
                }
            }
        }

        if(supportWidth > 0)
        {
            builder.anchoringSupportDistance((orientation == -1 ? -ciLeft : ciRight) + supportWidth);
        }
        else
        {
            builder.anchoringSupportDistance(0);
        }

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
            Genotype geno = context.getGenotype(referenceOrdinal);
            if(geno.hasExtendedAttribute(SV_FRAGMENT_COUNT) || geno.hasExtendedAttribute(
                    SGL_FRAGMENT_COUNT))
            {
                Integer var = asInteger(geno.getExtendedAttribute(context.hasAttribute(PAR_ID) | context.hasAttribute(MATE_ID)
                        ? SV_FRAGMENT_COUNT
                        : SGL_FRAGMENT_COUNT));
                Integer ref = asInteger(geno.getExtendedAttribute(REF_READ_COVERAGE));
                Integer refpair = asInteger(geno.getExtendedAttribute(REF_READPAIR_COVERAGE));
                builder.normalVariantFragmentCount(var);
                builder.normalReferenceFragmentCount(ref + (ignoreRefpair ? 0 : refpair));
            }
        }

        if(context.getGenotype(tumorOrdinal) != null)
        {
            Genotype geno = context.getGenotype(tumorOrdinal);
            if(geno.hasExtendedAttribute(SV_FRAGMENT_COUNT) || geno.hasExtendedAttribute(
                    SGL_FRAGMENT_COUNT))
            {
                Integer var = asInteger(geno.getExtendedAttribute(context.hasAttribute(PAR_ID) | context.hasAttribute(MATE_ID)
                        ? SV_FRAGMENT_COUNT
                        : SGL_FRAGMENT_COUNT));
                Integer ref = asInteger(geno.getExtendedAttribute(REF_READ_COVERAGE));
                Integer refpair = asInteger(geno.getExtendedAttribute(REF_READPAIR_COVERAGE));
                builder.tumorVariantFragmentCount(var);
                builder.tumorReferenceFragmentCount(ref + (ignoreRefpair ? 0 : refpair));
            }
        }

        return builder;
    }

    private static Double unadjustedAllelicFrequency(final VariantContext context)
    {
        if(context.hasAttribute(TAF))
            return context.getAttributeAsDoubleList(TAF, 0.0).get(0);

        return null;
    }

    private static Integer asInteger(Object obj)
    {
        if(obj == null)
        {
            return null;
        }
        else if(obj instanceof Integer)
        {
            return (Integer) obj;
        }
        else
        {
            final String strObj = obj.toString();
            if(strObj == null || strObj.isEmpty())
            {
                return null;
            }
            else
            {
                return Integer.parseInt(strObj);
            }
        }
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

    private static boolean imprecise(final VariantContext context)
    {
        final String impreciseStr = context.getAttributeAsString(IMPRECISE, "");
        boolean isPrecise = impreciseStr.isEmpty() || !impreciseStr.equals("true");
        return !isPrecise;
    }

    private static StructuralVariantType type(final VariantContext context)
    {
        return StructuralVariantType.fromAttribute((String) context.getAttribute(SVTYPE));
    }
}
