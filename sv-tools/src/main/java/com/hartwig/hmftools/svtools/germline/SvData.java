package com.hartwig.hmftools.svtools.germline;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INS;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.seIndex;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.BQ;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.BVF;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.HOMSEQ;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.QUAL;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.VF;

import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantType;

import htsjdk.variant.variantcontext.VariantContext;

public class SvData
{
    // VCF data
    private final String mId;
    private final StructuralVariantType mType;
    private final String[] mChromosome;
    private final int[] mPosition;
    private final byte[] mOrientation;

    private final int mReferenceOrdinal;
    private final int mTumorOrdinal;
    private VariantContext[] mContexts;

    private int mPonCount;
    private String mHotspotGeneInfo;

    // classification state
    private String[] mAsmbSvIds;

    public SvData(
            final String id, final StructuralVariantType type,
            final String chrStart, final String chrEnd, int posStart, int posEnd, byte orientStart, byte orientEnd,
            final VariantContext contextStart, final VariantContext contextEnd, final int referenceOrdinal, final int tumorOrdinal)
    {
        mId = id;
        mType = type;
        mChromosome = new String[] { chrStart, chrEnd };
        mPosition = new int[] { posStart, posEnd };
        mOrientation = new byte[] { orientStart, orientEnd };
        mReferenceOrdinal = referenceOrdinal;
        mTumorOrdinal = tumorOrdinal;
        mContexts = new VariantContext[] { contextStart, contextEnd };

        mPonCount = 0;
        mHotspotGeneInfo = "";
        mAsmbSvIds = new String[SE_PAIR];
    }

    public static SvData from(final StructuralVariant sv, final GenotypeIds genotypeIds)
    {
        boolean sglBreakend = sv.type() == SGL;

        return new SvData(
                sv.id(), sv.type(), sv.chromosome(true), !sglBreakend ? sv.chromosome(false) : "0",
                sv.position(true).intValue(), !sglBreakend ? sv.position(false).intValue() : -1,
                sv.orientation(true), !sglBreakend ? sv.orientation(false) : 0,
                sv.startContext(), sv.endContext(), genotypeIds.ReferenceOrdinal, genotypeIds.TumorOrdinal);
    }

    public String id() { return mId; }

    public String[] chromosomes() { return mChromosome; }
    public String chromosomeStart() { return mChromosome[SE_START]; }
    public String chromosomeEnd() { return mChromosome[SE_END]; }

    public int[] positions() { return mPosition; }
    public int posStart() { return mPosition[SE_START]; }
    public int posEnd() { return mPosition[SE_END]; }

    public byte[] orientations() { return mOrientation; }
    public byte orientStart() { return mOrientation[SE_START]; }
    public byte orientEnd() { return mOrientation[SE_END]; }

    public StructuralVariantType type() { return mType; }
    public boolean isSgl() { return mType == SGL; }

    public VariantContext[] contexts() { return mContexts; }
    public VariantContext contextStart() { return mContexts[SE_START]; }
    public VariantContext contextEnd() { return mContexts[SE_END]; }

    public double qual()
    {
        if(isSgl())
            return mContexts[SE_START].getAttributeAsDouble(BQ, 0);
        else
            return mContexts[SE_START].getAttributeAsDouble(QUAL, 0);
    }

    public int normalSupport()
    {
        if(isSgl())
            return mContexts[SE_START].getAttributeAsInt(BVF, 0);
        else
            return mContexts[SE_START].getAttributeAsInt(VF, 0);
    }

    public void setPonCount(int count) { mPonCount = count; }
    public int getPonCount() { return mPonCount; }

    public void setHotspotGeneInfo(final String info) { mHotspotGeneInfo = info; }
    public String getHotspotGeneInfo() { return mHotspotGeneInfo; }

    public boolean hasLength() { return hasLength(mType); }

    public static boolean hasLength(final StructuralVariantType type)
    {
        return type == INS || type == INV || type == DEL || type == DUP;
    }

    public int length() { return hasLength(mType) ? posEnd() - posStart() : 0; }

    public static int length(final StructuralVariant sv)
    {
        if(hasLength(sv.type()))
            return (int)(sv.position(false) - sv.position(true));

        return 0;
    }

    public String startHomology() { return contextStart().getAttributeAsString(HOMSEQ, ""); }
    public String endHomology() { return contextEnd() != null ? contextEnd().getAttributeAsString(HOMSEQ, "") : ""; }

    /*
        sampleName, sv.id(), sv.filter(), getDoubleValue(normalGenotype, QUAL), sv.type(),
        sv.chromosome(true), !sglBreakend ? sv.chromosome(false) : "0",
        sv.position(true), !sglBreakend ? sv.position(false) : -1,
        sv.orientation(true), !sglBreakend ? sv.orientation(false) : 0,
        getIntValue(normalGenotype, REF), getIntValue(normalGenotype, RP), getDoubleValue(normalGenotype, RPQ),
        getIntValue(normalGenotype, SR), getDoubleValue(normalGenotype, SRQ), getIntValue(normalGenotype, VF),
        sv.insertSequence(), homology, "", "");

                final String homology = variant.getAttributeAsString(HOMSEQ, "");


     */

    public String assemblySvIds(boolean isStart) { return mAsmbSvIds[seIndex(isStart)]; }
    public void setAssemblySvId(int seIndex, final String svId) { mAsmbSvIds[seIndex] = svId; }

    public String toString()
    {
        return String.format("%s:%s %s:%d - %s:%d",
                mId, mType.toString(), mChromosome[SE_START], mPosition[SE_START], mChromosome[SE_END], mPosition[SE_END]);
    }

    /*
    public final int NormalREF;
    public final int NormalRP;
    public final double NormalRPQ;
    public final int NormalSR;
    public final double NormalSRQ;
    public final int NormalVF;

     */

    /*
    public void createSV()
    {
        StructuralVariantLeg start = ImmutableStructuralVariantLegImpl.builder()
                .chromosome(ChrStart)
                .position(PosStart)
                .orientation((byte)OrientStart)
                .homology(Homology)
                .anchoringSupportDistance(0)
                .build();

        StructuralVariantLeg end = Type != SGL ? ImmutableStructuralVariantLegImpl.builder()
                .chromosome(ChrEnd)
                .position(PosEnd)
                .orientation((byte)OrientEnd)
                .homology(Homology)
                .anchoringSupportDistance(0)
                .build() : null;

        mSV = ImmutableStructuralVariantImpl.builder()
                .id(Id)
                .type(Type)
                .start(start)
                .end(end)
                .qualityScore(QualScore)
                .insertSequence(InsertSequence)
                .recovered(false)
                .hotspot(false)
                .build();
    }
    */

}
