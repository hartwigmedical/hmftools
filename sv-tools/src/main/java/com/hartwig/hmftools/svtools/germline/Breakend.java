package com.hartwig.hmftools.svtools.germline;

import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.BREAKEND_REGEX;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.CIPOS;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INS;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.svtools.germline.FilterConstants.SHORT_CALLING_SIZE;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.BQ;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.BVF;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.CIRPOS;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.QUAL;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.REF;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.REFPAIR;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.VF;

import java.util.List;
import java.util.regex.Matcher;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.sv.StructuralVariantLeg;
import com.hartwig.hmftools.common.sv.StructuralVariantType;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

public class Breakend
{
    public final String VcfId;
    public final VariantContext Context;
    public final String Chromosome;
    public final int Position;
    public final byte Orientation;

    public final double Qual;
    public final int TumorFragments;
    public final int ReferenceFragments;
    public final int ReferenceReads;
    public final int ReferencePairReads;

    public final Genotype RefGenotype;
    public final Genotype TumorGenotype;

    public final String Ref;
    public final String Alt;
    public final String InsertSequence;
    public final String OtherChromosome;
    public final int OtherPosition;
    public final byte OtherOrientation;

    public final String AssemblyInfo;

    public final Interval ConfidenceInterval;
    public final Interval RemoteConfidenceInterval;

    private final List<FilterType> mFilters;

    public Breakend(
            final String vcfId, final VariantContext context, final String chromosome, final int position, final byte orientation,
            final Genotype refGenotype, final Genotype tumorGenotype, final double qual, final int tumorFragments,
            final int refFrags, final int refReads, final int refPairReads, final String assemblyInfo)
    {
        VcfId = vcfId;
        Context = context;
        Chromosome = chromosome;
        Position = position;
        Orientation = orientation;

        RefGenotype = refGenotype;
        TumorGenotype = tumorGenotype;
        Qual = qual;
        TumorFragments = tumorFragments;
        ReferenceFragments = refFrags;
        ReferenceReads = refReads;
        ReferencePairReads = refPairReads;
        AssemblyInfo = assemblyInfo;

        ConfidenceInterval = VcfUtils.confidenceInterval(context, CIPOS);
        RemoteConfidenceInterval = VcfUtils.confidenceInterval(context, CIRPOS);

        final String[] refAltInsSeqParts = parseRefAlt(context);
        Ref = refAltInsSeqParts[REF_PART];
        Alt = refAltInsSeqParts[ALT_PART];
        InsertSequence = refAltInsSeqParts[INS_SEQ_PART];
        OtherChromosome = refAltInsSeqParts[CHR_PART];
        OtherPosition = Integer.parseInt(refAltInsSeqParts[POS_PART]);
        OtherOrientation = Byte.parseByte(refAltInsSeqParts[ORIENT_PART]);

        mFilters = Lists.newArrayList();
    }

    public static Breakend from(
            final String vcfId, final StructuralVariantType type, final StructuralVariantLeg svLeg, final VariantContext variantContext,
            final int referenceOrdinal, final int tumorOrdinal)
    {
        final Genotype tumorGenotype = variantContext.getGenotype(tumorOrdinal);
        final Genotype refGenotype = referenceOrdinal > 0 ? variantContext.getGenotype(referenceOrdinal) : null;

        final String qualTag = type == SGL ? BQ : QUAL;
        final String fragsTag = type == SGL ? BVF : VF;
        double qual = VcfUtils.getGenotypeAttributeAsDouble(tumorGenotype, qualTag, 0);

        int refFrags = 0;
        int refReads = 0;
        int refPairReads = 0;

        if(refGenotype != null)
        {
            refFrags = VcfUtils.getGenotypeAttributeAsInt(refGenotype, fragsTag, 0);
            refReads = VcfUtils.getGenotypeAttributeAsInt(refGenotype, REF, 0);
            refPairReads = VcfUtils.getGenotypeAttributeAsInt(refGenotype, REFPAIR, 0);
        }

        int tumorFrags = VcfUtils.getGenotypeAttributeAsInt(tumorGenotype, fragsTag, 0);

        return new Breakend(
                vcfId, variantContext, svLeg.chromosome(), (int)svLeg.position(), svLeg.orientation(), refGenotype, tumorGenotype,
                qual, tumorFrags, refFrags, refReads, refPairReads, "");
    }

    public static Breakend realigned(final Breakend original, final VariantContext newContext, final int newPosition)
    {
        return new Breakend(
                original.VcfId, newContext, original.Chromosome, newPosition, original.Orientation, original.RefGenotype,
                original.TumorGenotype, original.Qual, original.TumorFragments, original.ReferenceFragments, original.ReferenceReads,
                original.ReferencePairReads, original.AssemblyInfo);
    }

    public int insertSequenceLength() { return 0; } // TODO

    public int minPosition() { return Position + ConfidenceInterval.Start; }
    public int maxPosition() { return Position + ConfidenceInterval.End; }

    public List<FilterType> getFilters() { return mFilters; }

    public void addFilter(final FilterType filter)
    {
        if(!mFilters.contains(filter))
            mFilters.add(filter);
    }

    private static final int REF_PART = 0;
    private static final int ALT_PART = 1;
    private static final int INS_SEQ_PART = 2;
    private static final int CHR_PART = 3;
    private static final int POS_PART = 4;
    private static final int ORIENT_PART = 5;

    public static String[] parseRefAlt(final VariantContext variantContext)
    {
        String[] refAltParts = {"", "", "", "", "", ""};

        final String ref = variantContext.getAlleles().get(0).getDisplayString();
        refAltParts[REF_PART] = ref;

        final String alt = variantContext.getAlleles().get(1).getDisplayString();

        if(alt.startsWith("."))
        {
            refAltParts[ALT_PART] = alt.substring(alt.length() - 1);
            refAltParts[INS_SEQ_PART] = alt.substring(ref.length(), alt.length() - 1);
        }
        else if(alt.endsWith("."))
        {
            refAltParts[ALT_PART] = alt.substring(0, 1);
            refAltParts[INS_SEQ_PART] = alt.substring(1, alt.length() - ref.length());
        }
        else
        {
            final Matcher match = BREAKEND_REGEX.matcher(alt);

            if(!match.matches())
                return refAltParts;

            if(match.group(1).length() > 0)
            {
                String initialSequence = match.group(1).substring(1);
                refAltParts[INS_SEQ_PART] = initialSequence.substring(ref.length());
                refAltParts[ALT_PART] = alt.substring(0, 1);
            }
            else
            {
                String finalSequence = match.group(4).substring(0, match.group(4).length() - 1);
                refAltParts[INS_SEQ_PART] = finalSequence.substring(0, finalSequence.length() - ref.length());
                refAltParts[ALT_PART] = alt.substring(alt.length() - 1);
            }

            refAltParts[ORIENT_PART] = match.group(2).equals("]") ? "1" : "-1";

            String[] chrPos = match.group(3).split(":");
            refAltParts[CHR_PART] = chrPos[0];
            refAltParts[POS_PART] = chrPos[1];
        }

        return refAltParts;
    }

    public static String formPairedAltString(
            final String alt, final String insertSequence, final String chromosome, int position, byte orientStart, byte orientEnd)
    {
        if(orientStart == POS_ORIENT && orientEnd == NEG_ORIENT)
            return String.format("%s%s[%s:%d[", alt, insertSequence, chromosome, position);
        else if(orientStart == POS_ORIENT && orientEnd == POS_ORIENT)
            return String.format("%s%s]%s:%d]", alt, insertSequence, chromosome, position);
        else if(orientStart == NEG_ORIENT && orientEnd == NEG_ORIENT)
            return String.format("[%s:%d[%s%s", chromosome, position, insertSequence, alt);
        else
            return String.format("]%s:%d]%s%s", chromosome, position, insertSequence, alt);
    }

    public static String formSingleAltString(final String alt, final String insertSequence, byte orientation)
    {
        if(orientation == POS_ORIENT)
            return String.format("%s%s.", alt, insertSequence);
        else
            return String.format(".%s%s", insertSequence, alt);
    }


    /*
    val startBreakend: Breakend = Breakend(contig, start + confidenceInterval.first, start + confidenceInterval.second, orientation)
    val endBreakend: Breakend? = (variantType as? Paired)?.let { Breakend(it.otherChromosome, it.otherPosition + remoteConfidenceInterval.first, it.otherPosition + remoteConfidenceInterval.second, it.endOrientation) }
     */


}
