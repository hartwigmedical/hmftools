package com.hartwig.hmftools.common.fusion;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.fusion.KnownFusionCache.KF_LOGGER;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

public class KnownFusionData
{
    public final KnownFusionType Type;
    public final String FiveGene;
    public final String ThreeGene;
    public final String CancerTypes;
    public final String PubMedId;

    private boolean mValidData;

    // type-specific data:
    public boolean mHighImpactPromiscuous;

    private int[] mDownstreamDistance; // used for some known pair 5' genes and IG 3' genes

    // exon deletion
    private String mSpecificExonsTransName;

    private int[] mFiveExonRange;
    private int[] mThreeExonRange;

    // IG region
    private ChrBaseRegion mIgRegion;
    private byte mIgStrand;

    // 3' gene alternative mappings
    private final List<ChrBaseRegion> mThreeGeneAltRegions;

    private static final String FLD_TYPE = "Type";
    private static final String FLD_FIVE_GENE = "FiveGene";
    private static final String FLD_THREE_GENE = "ThreeGene";
    private static final String FLD_PUB_MED = "PubMedId";
    private static final String FLD_CANCER_TYPES = "CancerTypes";
    private static final String FLD_KNOWN_EXON_TRANS = "KnownExonTranscript";
    private static final String FLD_KNOWN_EXON_UP_RANGE = "KnownExonUpRange";
    private static final String FLD_KNOWN_EXON_DOWN_RANGE = "KnownExonDownRange";
    private static final String FLD_HIGH_IMPACT_PROM = "HighImpactPromiscuous";
    private static final String FLD_OVERRIDES = "Overrides";

    public static final String OVERRIDE_IG_RANGE = "IG_RANGE";
    public static final String OVERRIDE_ALTS = "ALTS";
    public static final String OVERRIDE_UP_DISTANCE = "UP_GENE_DOWNSTREAM_DISTANCE";
    public static final String OVERRIDE_DOWN_DISTANCE = "DOWN_GENE_DOWNSTREAM_DISTANCE";
    public static final String ALT_DATA = "ALT";

    private static final String FILE_DELIM = ",";
    private static final String ITEM_DELIM = ";";
    private static final String OVERRIDES_DELIM = " ";
    private static final String OVERRIDES_ID_DELIM = "=";

    public KnownFusionData(
            final KnownFusionType type, final String fiveGene, final String threeGene, final String cancerTypes, final String pubMedId)
    {
        Type = type;
        FiveGene = fiveGene;
        ThreeGene = threeGene;
        CancerTypes = cancerTypes;
        PubMedId = pubMedId;

        mHighImpactPromiscuous = false;
        mSpecificExonsTransName = "";
        mFiveExonRange = new int[SE_PAIR];
        mThreeExonRange = new int[SE_PAIR];
        mIgRegion = null;
        mIgStrand = 0;
        mDownstreamDistance = new int[] {0, 0};
        mThreeGeneAltRegions = Lists.newArrayList();

        mValidData = true;
    }

    public void setInvalid() { mValidData = false; }
    public boolean validData() { return mValidData; }

    public static KnownFusionData fromCsv(final String data, final Map<String,Integer> fieldIndexMap)
    {
        final String[] items = data.split(FILE_DELIM, -1);

        KnownFusionData kfData = new KnownFusionData(
                KnownFusionType.valueOf(items[fieldIndexMap.get(FLD_TYPE)]),
                items[fieldIndexMap.get(FLD_FIVE_GENE)],
                items[fieldIndexMap.get(FLD_THREE_GENE)],
                items[fieldIndexMap.get(FLD_CANCER_TYPES)],
                items[fieldIndexMap.get(FLD_PUB_MED)]);

        final String knownExonTrans = items[fieldIndexMap.get(FLD_KNOWN_EXON_TRANS)];
        final String knownExonUpRange = items[fieldIndexMap.get(FLD_KNOWN_EXON_UP_RANGE)];
        final String knownExonDownRange = items[fieldIndexMap.get(FLD_KNOWN_EXON_DOWN_RANGE)];

        try
        {
            if(!knownExonTrans.isEmpty())
            {
                kfData.setKnownExonData(knownExonTrans, knownExonUpRange, knownExonDownRange);
            }

            if(items[fieldIndexMap.get(FLD_HIGH_IMPACT_PROM)].equalsIgnoreCase("TRUE"))
                kfData.setHighImpactPromiscuous();

            final String overrides = items[fieldIndexMap.get(FLD_OVERRIDES)];

            if(!overrides.isEmpty())
                kfData.applyOverrides(overrides);
        }
        catch(Exception e)
        {
            KF_LOGGER.error("failed to parse specific data for known fusion({}): error({})", kfData, e.toString());
        }

        return kfData;
    }

    public void setKnownExonData(final String knownExonTrans, final String knownExonUpRange, final String knownExonDownRange)
    {
        mSpecificExonsTransName = knownExonTrans;

        if(!knownExonUpRange.isEmpty())
        {
            final String[] exons = knownExonUpRange.split(ITEM_DELIM);
            mFiveExonRange[SE_START]= Integer.parseInt(exons[SE_START]);
            mFiveExonRange[SE_END]= Integer.parseInt(exons[SE_END]);
        }

        if(!knownExonDownRange.isEmpty())
        {
            final String[] exons = knownExonDownRange.split(ITEM_DELIM);
            mThreeExonRange[SE_START]= Integer.parseInt(exons[SE_START]);
            mThreeExonRange[SE_END]= Integer.parseInt(exons[SE_END]);
        }
    }

    public void applyOverrides(final String overrides)
    {
        for(final String overrideItem : overrides.split(OVERRIDES_DELIM))
        {
            final String overrideName = overrideItem.split(OVERRIDES_ID_DELIM)[0];
            final String overrideData = overrideItem.split(OVERRIDES_ID_DELIM)[1];

            if(overrideName.equals(OVERRIDE_UP_DISTANCE))
            {
               mDownstreamDistance[FS_UP] = Integer.parseInt(overrideData);
            }
            else if(overrideName.equals(OVERRIDE_DOWN_DISTANCE))
            {
                mDownstreamDistance[FS_DOWN] = Integer.parseInt(overrideData);
            }
            else if(overrideName.equals(OVERRIDE_ALTS))
            {
                final String[] altItems = overrideData.split(ITEM_DELIM);
                int index = 0;
                while(index < altItems.length)
                {
                    if(!altItems[index].equals(ALT_DATA) || altItems.length - index < 4)
                    {
                        setInvalid();
                        return;
                    }

                    ++index;

                    mThreeGeneAltRegions.add(
                            new ChrBaseRegion(altItems[index], Integer.parseInt(altItems[index + 1]), Integer.parseInt(altItems[index + 2])));
                    index += 3;
                }
            }
            else if(overrideName.equals(OVERRIDE_IG_RANGE))
            {
                final String[] igRangeItems = overrideData.split(ITEM_DELIM);
                mIgStrand = Byte.parseByte(igRangeItems[0]);
                mIgRegion = new ChrBaseRegion(igRangeItems[1], Integer.parseInt(igRangeItems[2]), Integer.parseInt(igRangeItems[3]));
            }
        }
    }

    public boolean isHighImpactPromiscuous() { return mHighImpactPromiscuous; }
    public void setHighImpactPromiscuous() { mHighImpactPromiscuous = true; }

    public int downstreamDistance(int stream) { return mDownstreamDistance[stream]; }

    public String specificExonsTransName() { return mSpecificExonsTransName; }

    public int[] fiveGeneExonRange() { return mFiveExonRange; }
    public int[] threeGeneExonRange() { return mThreeExonRange; }

    public ChrBaseRegion igRegion() { return mIgRegion; }

    public boolean withinIgRegion(final String chromosome, int position)
    {
        return mIgRegion != null && mIgRegion.containsPosition(chromosome, position);
    }

    public boolean matchesIgGene(final String chromosome, int position, byte orientation)
    {
        return mIgStrand == orientation && withinIgRegion(chromosome, position);
    }

    public final List<ChrBaseRegion> getThreeGeneAltRegions() { return mThreeGeneAltRegions; }

    public String toString()
    {
        return String.format("%s: genes(%s - %s) ct(%s)",
                Type, FiveGene, ThreeGene, CancerTypes);
    }

    /*
    private void setTypeInfo()
    {
        if(Type == KNOWN_PAIR)
        {
            if(!mOtherData.isEmpty())
                mDownstreamDistance = Integer.parseInt(mOtherData);
        }
        else if(Type == EXON_DEL_DUP)
        {
            final String[] items = mOtherData.split(ITEM_DELIM);
            if(items.length != 5)
            {
                mValidData = false;
                return;
            }

            mSpecificExonsTransName = items[0];
            mMinFusedExons[FS_UPSTREAM] = Integer.parseInt(items[1]);
            mMaxFusedExons[FS_UPSTREAM] = Integer.parseInt(items[2]);
            mMinFusedExons[FS_DOWNSTREAM] = Integer.parseInt(items[3]);
            mMaxFusedExons[FS_DOWNSTREAM] = Integer.parseInt(items[4]);
        }
        else if(Type == IG_KNOWN_PAIR || Type == IG_PROMISCUOUS)
        {
            final String[] items = mOtherData.split(ITEM_DELIM);

            if(items.length != 5)
            {
                mValidData = false;
                return;
            }

            mIgStrand = Byte.parseByte(items[0]);

            mIgRegion = new SvRegion(items[1], Integer.parseInt(items[2]), Integer.parseInt(items[3]));
            mDownstreamDistance = Integer.parseInt(items[4]);
        }
        else if(Type == KNOWN_PAIR_UNMAPPABLE_3)
        {
            final String[] items = mOtherData.split(ITEM_DELIM);

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
                mDownstreamDistance = Integer.parseInt(items[4]);

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

    */

}
