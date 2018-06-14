package com.hartwig.hmftools.common.variant.structural;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import com.google.common.base.Preconditions;
import com.google.common.base.Strings;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.filter.ChromosomeFilter;

import htsjdk.samtools.util.StringUtil;
import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.CompoundFilter;
import htsjdk.variant.variantcontext.filter.PassingVariantFilter;
import htsjdk.variant.variantcontext.filter.VariantContextFilter;

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
                    String mate = (String) context.getAttribute(MATE_ID);
                    if (mate == null) {
                        mate = (String) context.getAttribute(PAR_ID);
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
        final String filtersStr = context.getFilters().toString();

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
        final List<Integer> cipos = context.getAttributeAsIntList(CIPOS, 0);
        final List<Integer> ihompos = context.getAttributeAsIntList(IHOMPOS, 0);
        final StructuralVariantLeg startLeg = ImmutableStructuralVariantLegImpl.builder()
                .chromosome(context.getContig())
                .position(start)
                .startPosition(start + (cipos.size() != 2 ? 0 : cipos.get(0)))
                .endPosition(start + (cipos.size() != 2 ? 0 : cipos.get(1)))
                .impreciseHomologyIntervalStart(start + (ihompos.size() != 2 ? 0 : ihompos.get(0)))
                .impreciseHomologyIntervalEnd(start + (ihompos.size() != 2 ? 0 : ihompos.get(1)))
                .orientation(startOrientation)
                .homology(context.getAttributeAsString(HOM_SEQ, ""))
                .alleleFrequency(af.size() == 2 ? af.get(0) : null)
                .build();

        final StructuralVariantLeg endLeg = ImmutableStructuralVariantLegImpl.builder()
                .chromosome(context.getContig())
                .position(end)
                .startPosition(end + (cipos.size() != 2 ? 0 : cipos.get(0)))
                .endPosition(end + (cipos.size() != 2 ? 0 : cipos.get(1)))
                .impreciseHomologyIntervalStart(end + (ihompos.size() != 2 ? 0 : ihompos.get(0)))
                .impreciseHomologyIntervalEnd(end + (ihompos.size() != 2 ? 0 : ihompos.get(1)))
                .orientation(endOrientation)
                .homology("")
                .alleleFrequency(af.size() == 2 ? af.get(1) : null)
                .build();

        return ImmutableStructuralVariantImpl.builder()
                .id(context.getID())
                .start(startLeg)
                .end(endLeg)
                .insertSequence(insertedSequence)
                .type(type)
                .filter(filtersStr)
                .imprecise(imprecise(context))
                .somaticScore(somaticScore)
                .qualityScore(context.getPhredScaledQual())
                .build();
    }

    @NotNull
    private static StructuralVariant create(@NotNull VariantContext first, @NotNull VariantContext second) {
        Preconditions.checkArgument(StructuralVariantType.BND.equals(type(first)));
        Preconditions.checkArgument(StructuralVariantType.BND.equals(type(second)));

        final int start = first.hasAttribute(BPI_START) ? first.getAttributeAsInt(BPI_START, -1) : first.getStart();
        final int end = second.hasAttribute(BPI_START) ? second.getAttributeAsInt(BPI_START, -1) : second.getStart();
        final List<Double> af = first.hasAttribute(BPI_AF) ? first.getAttributeAsDoubleList(BPI_AF, 0.0) : Collections.emptyList();
        final String filtersStr = first.getFilters().toString() + ";" + second.getFilters().toString();

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

        final int somaticScore = first.getAttributeAsInt(SOMATIC_SCORE, 0);

        final List<Integer> firstcipos = first.getAttributeAsIntList(CIPOS, 0);
        final List<Integer> secondcipos = second.getAttributeAsIntList(CIPOS, 0);
        final List<Integer> firstihompos = first.getAttributeAsIntList(IHOMPOS, 0);
        final List<Integer> secondihompos = first.getAttributeAsIntList(IHOMPOS, 0);

        final StructuralVariantLeg startLeg = ImmutableStructuralVariantLegImpl.builder()
                .chromosome(first.getContig())
                .position(start)
                .startPosition(start + (firstcipos.size() != 2 ? 0 : firstcipos.get(0)))
                .endPosition(start + (firstcipos.size() != 2 ? 0 : firstcipos.get(1)))
                .impreciseHomologyIntervalStart(start + (firstihompos.size() != 2 ? 0 : firstihompos.get(0)))
                .impreciseHomologyIntervalEnd(start + (firstihompos.size() != 2 ? 0 : firstihompos.get(1)))
                .orientation(startOrientation)
                .homology(first.getAttributeAsString(HOM_SEQ, ""))
                .alleleFrequency(af.size() == 2 ? af.get(0) : null)
                .build();

        final StructuralVariantLeg endLeg = ImmutableStructuralVariantLegImpl.builder()
                .chromosome(second.getContig())
                .position(end)
                .startPosition(end + (secondcipos.size() != 2 ? 0 : secondcipos.get(0)))
                .endPosition(end + (secondcipos.size() != 2 ? 0 : secondcipos.get(1)))
                .impreciseHomologyIntervalStart(end + (secondihompos.size() != 2 ? 0 : secondihompos.get(0)))
                .impreciseHomologyIntervalEnd(end + (secondihompos.size() != 2 ? 0 : secondihompos.get(1)))
                .orientation(endOrientation)
                .homology(second.getAttributeAsString(HOM_SEQ, ""))
                .alleleFrequency(af.size() == 2 ? af.get(1) : null)
                .build();

        return ImmutableStructuralVariantImpl.builder()
                .id(first.getID())
                .start(startLeg)
                .end(endLeg)
                .mateId(second.getID())
                .insertSequence(Strings.isNullOrEmpty(mantaInsertedSequence) ? insertedSequence : mantaInsertedSequence)
                .type(StructuralVariantType.BND)
                .filter(filtersStr)
                .imprecise(imprecise(first))
                .somaticScore(somaticScore)
                .qualityScore(first.getPhredScaledQual())
                .build();
    }
    @NotNull
    private static StructuralVariant createSingleBreakend(@NotNull VariantContext context) {
        Preconditions.checkArgument(StructuralVariantType.BND.equals(type(context)));
        Preconditions.checkArgument(singleBreakendRegex.matcher(context.getAlternateAllele(0).getDisplayString()).matches());

        final List<Double> af = context.hasAttribute(BPI_AF) ? context.getAttributeAsDoubleList(BPI_AF, 0.0) : Collections.emptyList();
        final String filtersStr = context.getFilters().toString();

        final String alt = context.getAlternateAllele(0).getDisplayString();
        // local orientation determined by the positioning of the anchoring bases
        final byte orientation = (byte)(alt.startsWith(".") ? -1 : 1);
        final int refLength = context.getReference().length();
        final String insertedSequence = orientation == -1 ? alt.substring(1, alt.length() - refLength) : alt.substring(refLength, alt.length() - 1);

        final List<Integer> cipos = context.getAttributeAsIntList(CIPOS, 0);

        final StructuralVariantLeg startLeg = ImmutableStructuralVariantLegImpl.builder()
                .chromosome(context.getContig())
                .position(context.getStart())
                .startPosition(context.getStart() + (cipos.size() != 2 ? 0 : cipos.get(0)))
                .endPosition(context.getStart() + (cipos.size() != 2 ? 0 : cipos.get(1)))
                .impreciseHomologyIntervalStart(context.getStart())
                .impreciseHomologyIntervalEnd(context.getStart())
                .orientation(orientation)
                .homology("")
                .alleleFrequency(af.size() >= 1 ? af.get(0) : null)
                .build();

        return ImmutableStructuralVariantImpl.builder()
                .id(context.getID())
                .start(startLeg)
                .insertSequence(insertedSequence)
                .type(StructuralVariantType.BND)
                .filter(filtersStr)
                .imprecise(imprecise(context))
                .somaticScore(context.getAttributeAsInt(SOMATIC_SCORE, 0))
                .qualityScore(context.getPhredScaledQual())
                .build();

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
