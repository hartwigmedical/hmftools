package com.hartwig.hmftools.neo.bind;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindCommon.AMINO_ACID_21ST;
import static com.hartwig.hmftools.neo.bind.BindScorer.INVALID_CALC;
import static com.hartwig.hmftools.neo.bind.RecognitionData.loadRecognitionData;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.aminoacid.BlosumMapping;

public class RecognitionSimilarity
{
    private final Map<String,Map<Integer,List<RecognitionData>>> mAllelePeptideDataMap; // counts data by peptide length>

    private final BlosumMapping mBlosumMapping;

    private boolean mCheckSelfSimilarity;
    private boolean mBlendAlleles;

    public RecognitionSimilarity()
    {
        mAllelePeptideDataMap = Maps.newHashMap();
        mBlosumMapping = new BlosumMapping();
        mCheckSelfSimilarity = false;
        mBlendAlleles = false;
    }

    public void setCheckSelfSimilarity(boolean toggle) { mCheckSelfSimilarity = toggle; }

    public boolean hasData() { return !mAllelePeptideDataMap.isEmpty(); }

    public double calcSimilarity(final String allele, final String peptide)
    {
        TopSimilarity topSimilarity = findTopSimilarity(allele, peptide);
        return topSimilarity != null ? topSimilarity.Similarity : INVALID_CALC;
    }

    private TopSimilarity findTopSimilarity(final String allele, final String peptide)
    {
        if(peptide.contains(AMINO_ACID_21ST))
            return null;

        Map<Integer,List<RecognitionData>> pepLenBindDataMap = mAllelePeptideDataMap.get(allele);

        if(pepLenBindDataMap == null)
            return null;

        List<RecognitionData> recogDataList = pepLenBindDataMap.get(peptide.length());

        if(recogDataList == null)
            return null;

        double topSimilarity = INVALID_CALC;
        String topPeptide = null;

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
                topPeptide = recogData.Peptide;

                if(topSimilarity == 0)
                    break;
            }
        }

        return topPeptide != null ? new TopSimilarity(topPeptide, topSimilarity) : null;
    }

    public double calcOtherAlleleSimilarity(final String allele, final String peptide)
    {
        double topSimilarity = INVALID_CALC;

        for(Map.Entry<String, Map<Integer, List<RecognitionData>>> alleleEntry : mAllelePeptideDataMap.entrySet())
        {
            String otherAllele = alleleEntry.getKey();

            if(otherAllele.equals(allele))
                continue;

            TopSimilarity topSim = findTopSimilarity(otherAllele, peptide);

            if(topSim != null)
            {
                if(topSimilarity == INVALID_CALC || topSim.Similarity < topSimilarity)
                    topSimilarity = topSim.Similarity;
            }
        }

        return topSimilarity;
    }

    public void logCrossAlleleSimilarities(final String filename)
    {
        boolean prevCheckSim = mCheckSelfSimilarity;
        mCheckSelfSimilarity = true;

        try
        {
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("Allele,Peptide,OtherAllele,TopPeptide,TopSimilarity");
            writer.newLine();

            for(Map.Entry<String, Map<Integer, List<RecognitionData>>> alleleEntry : mAllelePeptideDataMap.entrySet())
            {
                String allele = alleleEntry.getKey();

                for(Map.Entry<Integer, List<RecognitionData>> pepLenEntry : alleleEntry.getValue().entrySet())
                {
                    for(RecognitionData recogData : pepLenEntry.getValue())
                    {
                        for(String otherAllele : mAllelePeptideDataMap.keySet())
                        {
                            TopSimilarity topSimilarity = findTopSimilarity(otherAllele, recogData.Peptide);

                            if(topSimilarity == null)
                                continue;

                            writer.write(String.format("%s,%s,%s,%s,%.1f",
                                    allele, recogData.Peptide, otherAllele, topSimilarity.Peptide, topSimilarity.Similarity));
                            writer.newLine();
                        }
                    }
                }
            }

            writer.close();
        }
        catch (IOException e)
        {
            NE_LOGGER.error("failed to write recoognition similarities: {}", e.toString());
        }

        mCheckSelfSimilarity = prevCheckSim;
    }

    private class TopSimilarity
    {
        public final String Peptide;
        public final double Similarity;

        public TopSimilarity(final String peptide, final double similarity)
        {
            Peptide = peptide;
            Similarity = similarity;
        }
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
