package com.hartwig.hmftools.neo.bind;

import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.DELIMITER;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;

import java.io.BufferedReader;
import java.io.File;
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
    public final double Affinity;
    public final String Source;

    // optional fields
    private double mPredictedAffinity; // eg from MCF or another tool, only used for comparative purposes
    private double mAffinityPercentile;
    private double mPresentationPercentile;

    private double mScore;
    private double mRankPercentile;

    public static final String DELIM = ",";
    public static final String RANDOM_SOURCE = "Random";

    public BindData(final String allele, String peptide, double affinity, final String source)
    {
        Allele = allele;
        Peptide = peptide;
        Affinity = affinity;
        Source = source;

        mPredictedAffinity = -1;
        mPresentationPercentile = -1;
        mAffinityPercentile = -1;

        mScore = 0;
        mRankPercentile = -1;
    }

    public int peptideLength() { return Peptide.length(); }
    public boolean isRandom() { return Source.equals(RANDOM_SOURCE); }
    public boolean isTraining() { return !isRandom(); }

    public void setPredictionData(double predictedAffinity, double affinityPercentile, double presentationPercentile)
    {
        mPredictedAffinity = predictedAffinity;
        mAffinityPercentile = affinityPercentile;
        mPresentationPercentile = presentationPercentile;
    }

    public double predictedAffinity() { return mPredictedAffinity; }
    public double affinityPercentile() { return mAffinityPercentile; }
    public double presentationPercentile() { return mPresentationPercentile; }

    public void setScoreData(double score, double rankPerc)
    {
        mScore = score;
        mRankPercentile = rankPerc;
    }

    public double score() { return mScore; }
    public double rankPercentile() { return mRankPercentile; }

    public String toString()
    {
        return String.format("allele(%s) pep(%s) affinity(%.1f pred=%.1f) source(%s)",
                Allele, Peptide, Affinity, mPredictedAffinity, Source);
    }

    /*
    public static BindData fromCsv(
            final String data, int alleleIndex, int peptideIndex, int affinityIndex, int predictedIndex, int otherInfoIndex)
    {
        final String[] items = data.split(DELIM, -1);

        String allele = items[alleleIndex].replaceAll("HLA-", "");
        allele = allele.replaceAll(":", "").replaceAll("\\*", "");

        return new BindData(
                allele, items[peptideIndex], Double.parseDouble(items[affinityIndex]),
                Double.parseDouble(items[predictedIndex]), items[otherInfoIndex]);
    }

    public static BindData fromCsv(
            final String data, int alleleIndex, int peptideIndex, int predAffinityIndex, double maxAffinity)
    {
        final String[] items = data.split(DELIM, -1);
        double predictedAffinity = Double.parseDouble(items[predAffinityIndex]);
        return new BindData(items[alleleIndex], items[peptideIndex], maxAffinity, predictedAffinity, RANDOM_SOURCE);
    }
    */

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

            int alleleIndex = fieldsIndexMap.get("Allele");
            int peptideIndex = fieldsIndexMap.get("Peptide");
            Integer sourceIndex = fieldsIndexMap.get("Source");

            boolean isPredictionOnly = !fieldsIndexMap.containsKey("Affinity") && fieldsIndexMap.containsKey("PredictedAffinity");
            int affinityIndex = isPredictionOnly ? fieldsIndexMap.get("PredictedAffinity") : fieldsIndexMap.get("Affinity");

            // optional fields
            Integer predictedIndex = fieldsIndexMap.get("PredictedAffinity");
            Integer affinityPercIndex = fieldsIndexMap.get("AffinityPercentile");
            Integer presentationPercIndex = fieldsIndexMap.get("PresentationPercentile");

            // Allele,Peptide,PredictedAffinity,AffinityPercentile,PresentationPercentile,Source
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

                double affinity = !isPredictionOnly ? Double.parseDouble(items[affinityIndex]) : 50000;

                double predictedAffinity = predictedIndex != null ? Double.parseDouble(items[predictedIndex]) : -1;
                double affinityPerc = affinityPercIndex != null ? Double.parseDouble(items[affinityPercIndex]) : -1;
                double presentationPerc = presentationPercIndex != null ? Double.parseDouble(items[presentationPercIndex]) : -1;

                BindData bindData = new BindData(allele, peptide, affinity, source);

                bindData.setPredictionData(predictedAffinity, affinityPerc, presentationPerc);

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
