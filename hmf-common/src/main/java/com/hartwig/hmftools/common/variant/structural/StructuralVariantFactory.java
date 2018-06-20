package com.hartwig.hmftools.common.variant.structural;

import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import com.google.common.base.Preconditions;
import com.google.common.base.Strings;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.filter.ChromosomeFilter;

import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Genotype;
import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.CompoundFilter;
import htsjdk.variant.variantcontext.filter.PassingVariantFilter;
import htsjdk.variant.variantcontext.filter.VariantContextFilter;
import org.jetbrains.annotations.Nullable;

public class StructuralVariantFactory {

    private final static String TYPE = "SVTYPE";
    private final static String MATE_ID = "MATEID";
    private final static String PAR_ID = "PARID";
    private final static String INS_SEQ = "SVINSSEQ";
    private final static String LEFT_INS_SEQ = "LEFT_SVINSSEQ";
    private final static String RIGHT_INS_SEQ = "RIGHT_SVINSSEQ";
    private final static String HOM_SEQ = "HOMSEQ";
    private final static String BPI_START = "BPI_START";
    private final static String BPI_END = "BPI_END";
    private final static String BPI_AF = "BPI_AF";
    private final static String ALT = "ALT";
    private final static String IMPRECISE = "IMPRECISE";
    private final static String SOMATIC_SCORE = "SOMATICSCORE";
    private final static String IHOMPOS = "IHOMPOS";
    private final static String CIPOS = "CIPOS";
    private final static String VARIANT_FRAGMENT_BREAKEND_COVERAGE = "VF";
    private final static String REFERENCE_BREAKEND_READ_COVERAGE = "REF";
    private final static String REFERENCE_BREAKEND_READPAIR_COVERAGE = "REFPAIR";
    private final static String EVENT = "EVENT";
    private final static String LINKED_BY = "LINKED_BY";
    /**
     * Must match the small deldup threshold in scripts/gridss/gridss.config.R
     */
    private final static int SMALL_DELDUP_SIZE = 1000;
    private final static int NORMAL_GENOTYPE_ORDINAL = 0;
    private final static int TUMOUR_GENOTYPE_ORDINAL = 1;
    private final static Pattern breakendRegex = Pattern.compile("^(.*)([\\[\\]])(.+)[\\[\\]](.*)$");
    private final static Pattern singleBreakendRegex = Pattern.compile("^(([.].*)|(.*[.]))$");

    @NotNull
    private final Map<String, VariantContext> unmatched = Maps.newHashMap();
    @NotNull
    private final List<StructuralVariant> results = Lists.newArrayList();
    @NotNull
    private final VariantContextFilter filter;

    public StructuralVariantFactory(boolean filterPassesOnly) {
        final CompoundFilter filter = new CompoundFilter(true);

        if(filterPassesOnly) {
            filter.add(new PassingVariantFilter());
        }

        filter.add(new ChromosomeFilter());
        this.filter = filter;
    }

    public void addVariantContext(@NotNull VariantContext context) {
        if (filter.test(context)) {
            final StructuralVariantType type = type(context);
            if (type.equals(StructuralVariantType.BND)) {
                final boolean isSingleBreakend = singleBreakendRegex.matcher(context.getAlternateAllele(0).getDisplayString()).matches();
                if (isSingleBreakend) {
                    results.add(createSingleBreakend(context));
                } else {
                    String mate = (String) context.getAttribute(PAR_ID);
                    if (mate == null) {
                        mate = (String) context.getAttribute(MATE_ID);
                    }
                    if (unmatched.containsKey(mate)) {
                        results.add(create(unmatched.remove(mate), context));
                    } else {
                        unmatched.put(context.getID(), context);
                    }
                }
            } else {
                results.add(create(context));
            }
        }
    }

    @NotNull
    public List<StructuralVariant> results() {
        return results;
    }

    @NotNull
    private static StructuralVariant create(@NotNull VariantContext context) {
        final StructuralVariantType type = type(context);
        Preconditions.checkArgument(!StructuralVariantType.BND.equals(type));

        final int start = context.hasAttribute(BPI_START) ? context.getAttributeAsInt(BPI_START, -1) : context.getStart();
        final int end = context.hasAttribute(BPI_END) ? context.getAttributeAsInt(BPI_END, -1) : context.getEnd();
        final List<Double> af = context.hasAttribute(BPI_AF) ? context.getAttributeAsDoubleList(BPI_AF, 0.0) : Collections.emptyList();

        byte startOrientation = 0, endOrientation = 0;
        switch (type) {
            case INV:
                if (context.hasAttribute("INV3")) {
                    startOrientation = endOrientation = 1;
                } else if (context.hasAttribute("INV5")) {
                    startOrientation = endOrientation = -1;
                }
                break;
            case DEL:
                startOrientation = 1;
                endOrientation = -1;
                break;
            case INS:
                startOrientation = 1;
                endOrientation = -1;
                break;
            case DUP:
                startOrientation = -1;
                endOrientation = 1;
                break;
        }

        final int somaticScore = context.getAttributeAsInt(SOMATIC_SCORE, 0);

        String insertedSequence = "";

        if(type == StructuralVariantType.INS)
        {
            final String leftInsertSeq = context.getAttributeAsString(LEFT_INS_SEQ, "");
            final String rightInsertSeq = context.getAttributeAsString(RIGHT_INS_SEQ, "");
            if(!leftInsertSeq.isEmpty() && !rightInsertSeq.isEmpty())
            {
                insertedSequence = leftInsertSeq + "|" + rightInsertSeq;
            }
            else
            {
                List<Allele> alleles = context.getAlleles();
                if(alleles.size() > 1) {
                    insertedSequence = alleles.get(1).toString();

                    // remove the ref base from the start
                    insertedSequence = insertedSequence.substring(1, insertedSequence.length());
                }
            }
        }
        else
        {
            insertedSequence = context.getAttributeAsString(INS_SEQ, "");
        }
        final boolean isSmallDelDup = (end - start) <= SMALL_DELDUP_SIZE && (type == StructuralVariantType.DEL || type == StructuralVariantType.DUP);
        final StructuralVariantLeg startLeg = setLegCommon(ImmutableStructuralVariantLegImpl.builder(), context, isSmallDelDup)
                .chromosome(context.getContig())
                .position(start)
                .orientation(startOrientation)
                .homology(context.getAttributeAsString(HOM_SEQ, ""))
                .alleleFrequency(af.size() == 2 ? af.get(0) : null)
                .build();

        final StructuralVariantLeg endLeg = setLegCommon(ImmutableStructuralVariantLegImpl.builder(), context, isSmallDelDup)
                .chromosome(context.getContig())
                .position(end)
                .orientation(endOrientation)
                .homology("")
                .alleleFrequency(af.size() == 2 ? af.get(1) : null)
                .build();

        return setCommon(ImmutableStructuralVariantImpl.builder(), context)
                .start(startLeg)
                .end(endLeg)
                .insertSequence(insertedSequence)
                .type(type)
                .filter(filters(context, null))
                .build();
    }

    @NotNull
    private static StructuralVariant create(@NotNull VariantContext first, @NotNull VariantContext second) {
        Preconditions.checkArgument(StructuralVariantType.BND.equals(type(first)));
        Preconditions.checkArgument(StructuralVariantType.BND.equals(type(second)));

        final int start = first.hasAttribute(BPI_START) ? first.getAttributeAsInt(BPI_START, -1) : first.getStart();
        final int end = second.hasAttribute(BPI_START) ? second.getAttributeAsInt(BPI_START, -1) : second.getStart();
        final List<Double> af = first.hasAttribute(BPI_AF) ? first.getAttributeAsDoubleList(BPI_AF, 0.0) : Collections.emptyList();

        final String alt = first.getAlternateAllele(0).getDisplayString();
        final Matcher match = breakendRegex.matcher(alt);
        if (!match.matches()) {
            throw new IllegalArgumentException(String.format("ALT %s is not in breakend notation", alt));
        }
        // local orientation determined by the positioning of the anchoring bases
        final byte startOrientation = (byte)(match.group(1).length() > 0 ? 1 : -1);
        // other orientation determined by the direction of the brackets
        final byte endOrientation = (byte)(match.group(2).equals("]") ? 1 : -1);
        // grab the inserted sequence by removing 1 base from the reference anchoring bases
        String insertedSequence = match.group(1).length() > 0 ?
                match.group(1).substring(1) :
                match.group(4).substring(0, match.group(4).length() - 1);

        final String mantaInsertedSequence = first.getAttributeAsString(INS_SEQ, "");
        final boolean isSmallDelDup = first.getContig().equals(second.getContig()) &&
                Math.abs(first.getStart() - second.getStart()) <= SMALL_DELDUP_SIZE &&
                startOrientation != endOrientation;

        final StructuralVariantLeg startLeg = setLegCommon(ImmutableStructuralVariantLegImpl.builder(), first, isSmallDelDup)
                .position(start)
                .orientation(startOrientation)
                .homology(first.getAttributeAsString(HOM_SEQ, ""))
                .alleleFrequency(af.size() == 2 ? af.get(0) : null)
                .build();

        final StructuralVariantLeg endLeg = setLegCommon(ImmutableStructuralVariantLegImpl.builder(), second, isSmallDelDup)
                .position(end)
                .orientation(endOrientation)
                .homology(second.getAttributeAsString(HOM_SEQ, ""))
                .alleleFrequency(af.size() == 2 ? af.get(1) : null)
                .build();

        return setCommon(ImmutableStructuralVariantImpl.builder(), first)
                .start(startLeg)
                .end(endLeg)
                .mateId(second.getID())
                .insertSequence(Strings.isNullOrEmpty(mantaInsertedSequence) ? insertedSequence : mantaInsertedSequence)
                .type(StructuralVariantType.BND)
                .filter(filters(first, second))
                .build();
    }
    @NotNull
    private static StructuralVariant createSingleBreakend(@NotNull VariantContext context) {
        Preconditions.checkArgument(StructuralVariantType.BND.equals(type(context)));
        Preconditions.checkArgument(singleBreakendRegex.matcher(context.getAlternateAllele(0).getDisplayString()).matches());

        final List<Double> af = context.hasAttribute(BPI_AF) ? context.getAttributeAsDoubleList(BPI_AF, 0.0) : Collections.emptyList();

        final String alt = context.getAlternateAllele(0).getDisplayString();
        // local orientation determined by the positioning of the anchoring bases
        final byte orientation = (byte)(alt.startsWith(".") ? -1 : 1);
        final int refLength = context.getReference().length();
        final String insertedSequence = orientation == -1 ? alt.substring(1, alt.length() - refLength) : alt.substring(refLength, alt.length() - 1);

        final StructuralVariantLeg startLeg = setLegCommon(ImmutableStructuralVariantLegImpl.builder(), context, false)
                .orientation(orientation)
                .homology("")
                .alleleFrequency(af.size() >= 1 ? af.get(0) : null)
                .build();

        return setCommon(ImmutableStructuralVariantImpl.builder(), context)
                .start(startLeg)
                .insertSequence(insertedSequence)
                .type(StructuralVariantType.BND)
                .filter(filters(context, null))
                .build();

    }
    private static ImmutableStructuralVariantImpl.Builder setCommon(@NotNull ImmutableStructuralVariantImpl.Builder builder, @NotNull VariantContext context) {
        return builder
                .id(context.getID())
                .event(context.getAttributeAsString(EVENT, null))
                .linkedBy(context.getAttributeAsString(LINKED_BY, ""))
                .imprecise(imprecise(context))
                .somaticScore(context.hasAttribute(SOMATIC_SCORE) ? context.getAttributeAsInt(SOMATIC_SCORE, 0) : null)
                .qualityScore(context.getPhredScaledQual());
    }
    private static ImmutableStructuralVariantLegImpl.Builder setLegCommon(@NotNull ImmutableStructuralVariantLegImpl.Builder builder, @NotNull VariantContext context, boolean ignoreRefpair) {
        builder.chromosome(context.getContig());
        builder.position(context.getStart());
        if (context.hasAttribute(CIPOS)) {
            final List<Integer> cipos = context.getAttributeAsIntList(CIPOS, 0);
            if (cipos.size() == 2) {
                builder.startOffset(cipos.get(0));
                builder.endOffset(cipos.get(1));
            }
        }
        if (context.hasAttribute(IHOMPOS)) {
            final List<Integer> ihompos = context.getAttributeAsIntList(IHOMPOS, 0);
            if (ihompos.size() == 2) {
                builder.inexactHomologyOffsetStart(ihompos.get(0));
                builder.inexactHomologyOffsetEnd(ihompos.get(1));
            }
        }
        if (context.getGenotype(NORMAL_GENOTYPE_ORDINAL) != null) {
            Genotype geno = context.getGenotype(NORMAL_GENOTYPE_ORDINAL);
            if (geno.hasExtendedAttribute(VARIANT_FRAGMENT_BREAKEND_COVERAGE)) {
                Integer var = asInteger(geno.getExtendedAttribute(VARIANT_FRAGMENT_BREAKEND_COVERAGE));
                Integer ref = asInteger(geno.getExtendedAttribute(REFERENCE_BREAKEND_READ_COVERAGE));
                Integer refpair = asInteger(geno.getExtendedAttribute(REFERENCE_BREAKEND_READPAIR_COVERAGE));
                builder = builder.normalVariantFragmentCount(var);
                builder = builder.normalReferenceFragmentCount(ref + (ignoreRefpair ? 0 : refpair));
            }
        }
        if (context.getGenotype(TUMOUR_GENOTYPE_ORDINAL) != null) {
            Genotype geno = context.getGenotype(TUMOUR_GENOTYPE_ORDINAL);
            if (geno.hasExtendedAttribute(VARIANT_FRAGMENT_BREAKEND_COVERAGE)) {
                Integer var = asInteger(geno.getExtendedAttribute(VARIANT_FRAGMENT_BREAKEND_COVERAGE));
                Integer ref = asInteger(geno.getExtendedAttribute(REFERENCE_BREAKEND_READ_COVERAGE));
                Integer refpair = asInteger(geno.getExtendedAttribute(REFERENCE_BREAKEND_READPAIR_COVERAGE));
                builder = builder.tumourVariantFragmentCount(var);
                builder = builder.tumourReferenceFragmentCount(ref + (ignoreRefpair ? 0 : refpair));
            }
        }
        return builder;
    }
    private static Integer asInteger(Object obj) {
        if (obj == null) {
            return null;
        } else if (obj instanceof Integer) {
            return (Integer)obj;
        } else {
            final String strObj = obj.toString();
            if (strObj == null || strObj.isEmpty()) {
                return null;
            } else {
                return Integer.parseInt(strObj);
            }
        }
    }
    @NotNull
    private static String filters(@NotNull VariantContext context, @Nullable VariantContext pairedContext) {
        final HashSet<String> filters = new HashSet<>(context.getFilters());
        if (pairedContext != null) {
            filters.addAll(pairedContext.getFilters());
        }
        if (filters.size() > 1) {
            // Doesn't pass if a filter is applied to either of the two records
            filters.remove("PASS");
        }
        // TODO Collectors string concatenation
        final String filtersStr = filters.stream().collect(Collectors.joining(";"));
        return filtersStr;
    }
    @NotNull
    private static boolean imprecise(@NotNull VariantContext context) {
        final String impreciseStr = context.getAttributeAsString(IMPRECISE, "");
        boolean isPrecise = impreciseStr.isEmpty() || !impreciseStr.equals("true");
        return !isPrecise;
    }
    @NotNull
    private static StructuralVariantType type(@NotNull VariantContext context) {
        return StructuralVariantType.fromAttribute((String) context.getAttribute(TYPE));
    }
}
