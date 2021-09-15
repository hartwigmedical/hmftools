package com.hartwig.hmftools.neo.bind;

import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindCommon.AMINO_ACID_21ST;
import static com.hartwig.hmftools.neo.bind.BindScorer.INVALID_CALC;
import static com.hartwig.hmftools.neo.bind.RecognitionData.loadRecognitionData;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class RecognitionSimilarity
{
    private final Map<String,Map<Integer,List<RecognitionData>>> mAllelePeptideDataMap; // counts data by peptide length>

    private final BlosumMapping mBlosumMapping;

    private boolean mCheckSelfSimilarity;

    public RecognitionSimilarity()
    {
        mAllelePeptideDataMap = Maps.newHashMap();
        mBlosumMapping = new BlosumMapping();
        mCheckSelfSimilarity = false;
    }

    public void setCheckSelfSimilarity(boolean toggle) { mCheckSelfSimilarity = toggle; }

    public boolean hasData() { return !mAllelePeptideDataMap.isEmpty(); }

    public double calcSimilarity(final String allele, final String peptide)
    {
        if(peptide.contains(AMINO_ACID_21ST))
            return INVALID_CALC;

        Map<Integer,List<RecognitionData>> pepLenBindDataMap = mAllelePeptideDataMap.get(allele);

        if(pepLenBindDataMap == null)
            return INVALID_CALC;

        List<RecognitionData> recogDataList = pepLenBindDataMap.get(peptide.length());

        if(recogDataList == null)
            return INVALID_CALC;

        double topSimilarity = INVALID_CALC;

        for(RecognitionData recogData : recogDataList)
        {
            if(mCheckSelfSimilarity && recogData.Peptide.equals(peptide) && recogData.Allele.equals(allele))
                continue;

            double similarity = 0;
            boolean skip = false;

            for(int i = 0; i < peptide.length(); ++i)
            {
                char aa1 = peptide.charAt(i);
                char aa2 = recogData.Peptide.charAt(i);

                try
                {
                    int bs1 = mBlosumMapping.selfMapping(aa1);
                    int bs2 = mBlosumMapping.selfMapping(aa2);
                    int map = mBlosumMapping.map(aa1, aa2);

                    similarity += (bs1 + bs2) * 0.5 - map;

                    if(topSimilarity >= 0 && similarity >= topSimilarity) // early exit if cannot be better
                    {
                        skip = true;
                        break;
                    }
                }
                catch(Exception e)
                {
                    NE_LOGGER.error("invalid index({}) of peptide({}) or recogPeptide({}) for recognition similarity calc",
                            i, peptide, recogData.Peptide);
                }
            }

            if(skip)
                continue;

            if(topSimilarity < 0 || similarity < topSimilarity)
            {
                topSimilarity = similarity;

                if(topSimilarity == 0)
                    break;
            }
        }

        return topSimilarity;
    }

    public boolean loadData(final String filename)
    {
        List<RecognitionData> recognitionData = Lists.newArrayList();

        if(!loadRecognitionData(filename, recognitionData))
            return false;

        for(RecognitionData recogData : recognitionData)
        {
            if(!recogData.Immunogenic)
                continue;

            Map<Integer,List<RecognitionData>> pepLenBindDataMap = mAllelePeptideDataMap.get(recogData.Allele);

            if(pepLenBindDataMap == null)
            {
                pepLenBindDataMap = Maps.newHashMap();
                mAllelePeptideDataMap.put(recogData.Allele, pepLenBindDataMap);
            }

            List<RecognitionData> recogDataList = pepLenBindDataMap.get(recogData.peptideLength());

            if(recogDataList == null)
            {
                recogDataList = Lists.newArrayList();
                pepLenBindDataMap.put(recogData.peptideLength(), recogDataList);
            }

            recogDataList.add(recogData);
        }

        return true;
    }

}
