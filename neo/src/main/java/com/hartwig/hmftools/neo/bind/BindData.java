package com.hartwig.hmftools.neo.bind;

import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.DELIMITER;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_AFFINITY;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_ALLELE;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_DOWN_FLANK;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_PEPTIDE;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_PRED_AFFINITY;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_PRES_SCORE;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_SOURCE;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_UP_FLANK;
import static com.hartwig.hmftools.neo.bind.BindCommon.RANDOM_SOURCE;
import static com.hartwig.hmftools.neo.bind.BindConstants.INVALID_SCORE;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class BindData
{
    public final String Allele;
    public final String Peptide;
    public final String UpFlank;
    public final String DownFlank;
    public final String Source;

    // optional fields
    private boolean mHasMeasuredAffinity;
    private double mMeasuredAffinity;

    private boolean mHasPredictionData;
    private double mPredictedAffinity; // eg from MCF or another tool, only used for comparative purposes
    private double mAffinityPercentile;
    private double mPresentationScore;
    private double mPresentationPercentile;

    private double mScore;
    private double mRankPercentile;
    private double mLikelihood;
    private double mLikelihoodRank;

    private double mGlobalScore;
    private double mGlobalRankPercentile;

    public BindData(final String allele, final String peptide, final String source)
    {
        this(allele, peptide, source, "", "");
    }

    public BindData(final String allele, final String peptide, final String source, final String upFlank, final String downFlank)
    {
        Allele = allele;
        Peptide = peptide;
        Source = source;
        UpFlank = upFlank;
        DownFlank = downFlank;

        mHasMeasuredAffinity = false;
        mMeasuredAffinity = -1;

        mHasPredictionData = false;
        mPredictedAffinity = -1;
        mPresentationScore = -1;
        mPresentationPercentile = -1;
        mAffinityPercentile = -1;

        mScore = INVALID_SCORE;
        mRankPercentile = INVALID_SCORE;
        mLikelihood = INVALID_SCORE;
        mLikelihoodRank = INVALID_SCORE;

        mGlobalScore = 0;
        mGlobalRankPercentile = -1;
    }

    public int peptideLength() { return Peptide.length(); }
    public boolean isRandom() { return Source.equals(RANDOM_SOURCE); }
    public boolean isTraining() { return !isRandom(); }

    public boolean hasFlanks() { return !UpFlank.isEmpty() || !DownFlank.isEmpty(); }

    public void setMeasuredAffinity(double affinity)
    {
        mHasMeasuredAffinity = true;
        mMeasuredAffinity = affinity;
    }

    public double measuredAffinity() { return mMeasuredAffinity; }
    public boolean hasMeasuredAffinity() { return mHasMeasuredAffinity; }

    public void setPredictionData(
            double predictedAffinity, double affinityPercentile, double presentationScore, double presentationPercentile)
    {
        mHasPredictionData = true;
        mPredictedAffinity = predictedAffinity;
        mAffinityPercentile = affinityPercentile;
        mPresentationScore = presentationScore;
        mPresentationPercentile = presentationPercentile;
    }

    public boolean hasPredictionData() { return mHasPredictionData; }
    public double predictedAffinity() { return mPredictedAffinity; }
    public double affinityPercentile() { return mAffinityPercentile; }
    public double presentationScore() { return mPresentationScore; }
    public double presentationPercentile() { return mPresentationPercentile; }

    public void setScoreData(double score, double rankPerc, double likelihood, double likelihoodRank)
    {
        mScore = score;
        mRankPercentile = rankPerc;
        mLikelihood = likelihood;
        mLikelihoodRank = likelihoodRank;
    }

    public double score() { return mScore; }
    public double rankPercentile() { return mRankPercentile; }
    public double likelihood() { return mLikelihood; }
    public double likelihoodRank() { return mLikelihoodRank; }

    public void setGlobalScoreData(double score, double rankPerc)
    {
        mGlobalScore = score;
        mGlobalRankPercentile = rankPerc;
    }

    public double globalScore() { return mGlobalScore; }
    public double globalRankPercentile() { return mGlobalRankPercentile; }

    public String toString()
    {
        return String.format("allele(%s) pep(%s) affinity(%.1f pred=%.1f) source(%s)",
                Allele, Peptide, mMeasuredAffinity, mPredictedAffinity, Source);
    }

    public static String cleanAllele(final String allele)
    {
        return allele.replaceAll("HLA-", "").replaceAll(":", "").replaceAll("\\*", "");
    }

    public static boolean loadBindData(
            final String filename, boolean expectExists, final List<String> restrictedAlleles, final List<Integer> restrictedLengths,
            final Map<String,Map<Integer,List<BindData>>> allelePeptideMap)
    {
        if(filename == null || !Files.exists(Paths.get(filename)))
            return !expectExists;

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));
            String header = fileReader.readLine();

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, DELIMITER);

            int alleleIndex = fieldsIndexMap.get(FLD_ALLELE);
            int peptideIndex = fieldsIndexMap.get(FLD_PEPTIDE);

            // all other fields are optional
            Integer upFlankIndex = fieldsIndexMap.get(FLD_UP_FLANK);
            Integer downFlankIndex = fieldsIndexMap.get(FLD_DOWN_FLANK);
            Integer sourceIndex = fieldsIndexMap.get(FLD_SOURCE);
            Integer affinityIndex = fieldsIndexMap.get(FLD_AFFINITY);
            Integer predictedIndex = fieldsIndexMap.get(FLD_PRED_AFFINITY);
            Integer affinityPercIndex = fieldsIndexMap.get("AffinityPercentile");
            Integer presentationScoreIndex = fieldsIndexMap.get(FLD_PRES_SCORE);
            Integer presentationPercIndex = fieldsIndexMap.get("PresentationPercentile");

            int alleleCount = 0;
            int itemCount = 0;

            String line = "";

            while((line = fileReader.readLine()) != null)
            {
                final String[] items = line.split(DELIMITER, -1);

                String allele = cleanAllele(items[alleleIndex]);

                if(!restrictedAlleles.isEmpty() && !restrictedAlleles.contains(allele))
                {
                    if(restrictedAlleles.size() == alleleCount)
                        break;

                    continue;
                }

                String peptide = items[peptideIndex];

                if(!restrictedLengths.isEmpty() && !restrictedLengths.contains(peptide.length()))
                    continue;

                if(peptide.contains("X"))
                    continue;

                ++itemCount;

                String source = sourceIndex != null ? items[sourceIndex] : RANDOM_SOURCE;
                String upFlank = upFlankIndex != null ? items[upFlankIndex] : "";
                String downFlank = downFlankIndex != null ? items[downFlankIndex] : "";

                BindData bindData = new BindData(allele, peptide, source, upFlank, downFlank);

                if(affinityIndex != null)
                {
                    bindData.setMeasuredAffinity(Double.parseDouble(items[affinityIndex]));
                }

                if(predictedIndex != null || affinityPercIndex != null || presentationScoreIndex != null)
                {
                    double predictedAffinity = predictedIndex != null ? Double.parseDouble(items[predictedIndex]) : -1;
                    double affinityPerc = affinityPercIndex != null ? Double.parseDouble(items[affinityPercIndex]) : -1;
                    double presentationScore = presentationScoreIndex != null ? Double.parseDouble(items[presentationScoreIndex]) : -1;
                    double presentationPerc = presentationPercIndex != null ? Double.parseDouble(items[presentationPercIndex]) : -1;
                    bindData.setPredictionData(predictedAffinity, affinityPerc, presentationScore, presentationPerc);
                }

                Map<Integer,List<BindData>> pepLenBindDataMap = allelePeptideMap.get(allele);

                if(pepLenBindDataMap == null)
                {
                    ++alleleCount;
                    pepLenBindDataMap = Maps.newHashMap();
                    allelePeptideMap.put(allele, pepLenBindDataMap);
                }

                List<BindData> bindDataList = pepLenBindDataMap.get(bindData.peptideLength());

                if(bindDataList == null)
                {
                    bindDataList = Lists.newArrayList();
                    pepLenBindDataMap.put(bindData.peptideLength(), bindDataList);
                }

                bindDataList.add(bindData);
            }

            NE_LOGGER.info("loaded {} alleles with {} bind data items from file({})",
                    allelePeptideMap.size(), itemCount, filename);
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to read binding data file: {}", e.toString());
            return false;
        }

        return true;
    }
}
