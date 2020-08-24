package com.hartwig.hmftools.common.fusion;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWNSTREAM;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_PAIR;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UPSTREAM;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.EXON_DEL_DUP;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.IG_KNOWN_PAIR;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.IG_PROMISCUOUS;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.KNOWN_PAIR_UNMAPPABLE_3;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.sv.SvRegion;

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
    private SvRegion mIgRegion;
    private byte mIgStrand;
    private int mIgDownstreamDistance;

    // 3' gene alternative mappings
    private final List<SvRegion> mThreeGeneAltRegions;

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
    private static final String ALT_DATA = "ALT";

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
        mIgRegion = null;
        mIgStrand = 0;
        mIgDownstreamDistance = 0;
        mThreeGeneAltRegions = Lists.newArrayList();

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

            mIgRegion = new SvRegion(items[1], Integer.parseInt(items[2]), Integer.parseInt(items[3]));
            mIgDownstreamDistance = Integer.parseInt(items[4]);
        }
        else if(Type == KNOWN_PAIR_UNMAPPABLE_3)
        {
            final String[] items = mOtherData.split(OTHER_DATA_DELIMITER);

            // non-IG example: ALT;GL0002281;20000;125000;ALT;4;190930000;191030000;ALT;10;135420000;135520000
            // IG example: 1;14;106032614;107288051;0;ALT;GL0002281;20000;125000;ALT;4;190930000;191030000;ALT;10;135420000;135520000

            if(items.length < 4)
            {
                mValidData = false;
                return;
            }

            int index = 0;
            if(!items[0].equals(ALT_DATA))
            {
                if(items.length < 5)
                {
                    mValidData = false;
                    return;
                }

                mIgStrand = Byte.parseByte(items[0]);

                mIgRegion = new SvRegion(items[1], Integer.parseInt(items[2]), Integer.parseInt(items[3]));
                mIgDownstreamDistance = Integer.parseInt(items[4]);

                index += 5;
            }

            while(index < items.length)
            {
                if(!items[index].equals(ALT_DATA) || items.length - index < 4)
                {
                    mValidData = false;
                    return;
                }

                ++index;

                mThreeGeneAltRegions.add(new SvRegion(items[index], Integer.parseInt(items[index+1]), Integer.parseInt(items[index+2])));
                index += 3;
            }
        }
    }

    public String specificTransName() { return mSpecificTransName; }
    public int[] minFusedExons() { return mMinFusedExons; }
    public int[] maxFusedExons() { return mMaxFusedExons; }
    public SvRegion igRegion() { return mIgRegion; }

    public boolean withinIgRegion(final String chromosome, int position)
    {
        return mIgRegion != null && mIgRegion.containsPosition(chromosome, position);
    }

    public boolean matchesIgGene(final String chromosome, int position, byte orientation)
    {
        return mIgStrand == orientation && withinIgRegion(chromosome, position);
    }

    public int igDownstreamDistance() { return mIgDownstreamDistance; }

    public final List<SvRegion> getThreeGeneAltRegions() { return mThreeGeneAltRegions; }

    public String toString()
    {
        return String.format("%s: genes(%s - %s) ct(%s) otherData(%s)",
                Type, FiveGene, ThreeGene, CancerTypes, mOtherData);
    }
}
