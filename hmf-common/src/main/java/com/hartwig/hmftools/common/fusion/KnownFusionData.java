package com.hartwig.hmftools.common.fusion;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWNSTREAM;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_PAIR;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UPSTREAM;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.EXON_DEL_DUP;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.IG_KNOWN_PAIR;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.IG_PROMISCUOUS;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

import java.util.Map;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;

public class KnownFusionData
{
    public final KnownFusionType Type;
    public final String FiveGene;
    public final String ThreeGene;
    public final String CancerTypes;
    public final String PubMedId;

    private boolean mValidData;
    private final String mOtherData;

    // exon deletion
    private String mSpecificTransName;
    private int[] mMinFusedExons;
    private int[] mMaxFusedExons;

    // IG region
    private int[] mIgRegion;
    private String mIgChromosome;
    private byte mIgStrand;
    private int mIgDownstreamDistance;

    public static final int FIVE_GENE = 0;
    public static final int THREE_GENE = 1;

    private static final String FILE_DELIMITER = ",";
    private static final String OTHER_DATA_DELIMITER = ";";
    private static final String FLD_TYPE = "Type";
    private static final String FLD_FIVE_GENE = "FiveGene";
    private static final String FLD_THREE_GENE = "ThreeGene";
    private static final String FLD_PUB_MED = "PubMedId";
    private static final String FLD_CANCER_TYPES = "CancerTypes";
    private static final String FLD_OTHER_DATA = "OtherData";

    public KnownFusionData(
            final KnownFusionType type, final String fiveGene, final String threeGene, final String cancerTypes,
            final String pubMedId, final String otherData)
    {
        Type = type;
        FiveGene = fiveGene;
        ThreeGene = threeGene;
        CancerTypes = cancerTypes;
        PubMedId = pubMedId;
        mOtherData = otherData;

        mSpecificTransName = "";
        mMinFusedExons = new int[FS_PAIR];
        mMaxFusedExons = new int[FS_PAIR];
        mIgRegion = new int[SE_PAIR];
        mIgChromosome = "";
        mIgStrand = 0;
        mIgDownstreamDistance = 0;

        mValidData = true;

        setTypeInfo();
    }

    public boolean validData() { return mValidData; }

    public static KnownFusionData fromCsv(final String data, final Map<String,Integer> fieldIndexMap)
    {
        final String[] items = data.split(FILE_DELIMITER, -1);

        return new KnownFusionData(
                KnownFusionType.valueOf(items[fieldIndexMap.get(FLD_TYPE)]),
                items[fieldIndexMap.get(FLD_FIVE_GENE)],
                items[fieldIndexMap.get(FLD_THREE_GENE)],
                items[fieldIndexMap.get(FLD_CANCER_TYPES)],
                items[fieldIndexMap.get(FLD_PUB_MED)],
                items[fieldIndexMap.get(FLD_OTHER_DATA)]);
    }

    private void setTypeInfo()
    {
        if(Type == EXON_DEL_DUP)
        {
            final String[] items = mOtherData.split(OTHER_DATA_DELIMITER);
            if(items.length != 5)
            {
                mValidData = false;
                return;
            }

            mSpecificTransName = items[0];
            mMinFusedExons[FS_UPSTREAM] = Integer.parseInt(items[1]);
            mMaxFusedExons[FS_UPSTREAM] = Integer.parseInt(items[2]);
            mMinFusedExons[FS_DOWNSTREAM] = Integer.parseInt(items[3]);
            mMaxFusedExons[FS_DOWNSTREAM] = Integer.parseInt(items[4]);
        }
        else if(Type == IG_KNOWN_PAIR || Type == IG_PROMISCUOUS)
        {
            final String[] items = mOtherData.split(OTHER_DATA_DELIMITER);
            if(items.length != 5)
            {
                mValidData = false;
                return;
            }

            mIgStrand = Byte.parseByte(items[0]);
            mIgChromosome = items[1];
            mIgRegion[SE_START] = Integer.parseInt(items[2]);
            mIgRegion[SE_END] = Integer.parseInt(items[3]);
            mIgDownstreamDistance = Integer.parseInt(items[4]);
        }
    }

    public String specificTransName() { return mSpecificTransName; }
    public int[] minFusedExons() { return mMinFusedExons; }
    public int[] maxFusedExons() { return mMaxFusedExons; }
    public int[] igRegion() { return mIgRegion; }

    public boolean withinIgRegion(final String chromosome, int position)
    {
        return mIgChromosome.equals(chromosome) && position >= mIgRegion[SE_START] && position <= mIgRegion[SE_END];
    }

    public boolean matchesIgGene(final String chromosome, int position, byte orientation)
    {
        return mIgStrand == orientation && withinIgRegion(chromosome, position);
    }

    public int igDownstreamDistance() { return mIgDownstreamDistance; }

    public String toString()
    {
        return String.format("%s: genes(%s - %s) ct(%s) otherData(%s)",
                Type, FiveGene, ThreeGene, CancerTypes, mOtherData);
    }
}
