package com.hartwig.hmftools.neo.bind;

import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.DELIMITER;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindCommon.AMINO_ACID_21ST;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_ALLELE;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_DOWN_FLANK;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_PEPTIDE;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_SOURCE;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_TPM;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_UP_FLANK;
import static com.hartwig.hmftools.neo.bind.BindCommon.cleanAllele;
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
    private final List<String> mOtherData;

    private Double mTPM;

    // scoring fields

    private double mScore;
    private double mFlankScore;
    private double mRankPercentile;
    private double mLikelihood;
    private double mExpLikelihood;
    private double mLikelihoodRank;

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

        mOtherData = Lists.newArrayList();
        mTPM = 0.0;

        mScore = INVALID_SCORE;
        mFlankScore = INVALID_SCORE;
        mRankPercentile = INVALID_SCORE;
        mLikelihood = INVALID_SCORE;
        mExpLikelihood = INVALID_SCORE;
        mLikelihoodRank = INVALID_SCORE;
    }

    public int peptideLength() { return Peptide.length(); }

    public boolean hasFlanks() { return !UpFlank.isEmpty() || !DownFlank.isEmpty(); }

    public final List<String> getOtherData() { return mOtherData; }

    public void setTPM(double tpm) { mTPM = tpm; }
    public boolean hasTPM() { return mTPM != null; }
    public double tpm() { return mTPM != null ? mTPM : 0; }

    public void setScoreData(
            double score, double flankScore, double rankPerc, double likelihood, double expLikelihood, double likelihoodRank)
    {
        mScore = score;
        mFlankScore = flankScore;
        mRankPercentile = rankPerc;
        mLikelihood = likelihood;
        mExpLikelihood = expLikelihood;
        mLikelihoodRank = likelihoodRank;
    }

    public double score() { return mScore; }
    public double flankScore() { return mFlankScore; }
    public double rankPercentile() { return mRankPercentile; }
    public double likelihood() { return mLikelihood; }
    public double expressionLikelihood() { return mExpLikelihood; }
    public double likelihoodRank() { return mLikelihoodRank; }

    public String toString()
    {
        return String.format("allele(%s) peptide(%s) source(%s)", Allele, Peptide, Source);
    }

    public static boolean loadBindData(
            final String filename, final List<Integer> restrictedLengths, final Map<String,Map<Integer,List<BindData>>> allelePeptideMap)
    {
        return loadBindData(filename, restrictedLengths, allelePeptideMap, Maps.newHashMap());
    }

    public static boolean loadBindData(
            final String filename, final List<Integer> restrictedLengths, final Map<String,Map<Integer,List<BindData>>> allelePeptideMap,
            final Map<String,Integer> otherColumns)
    {
        if(filename == null || !Files.exists(Paths.get(filename)))
        {
            NE_LOGGER.error("binding data file({}) not found", filename);
            return false;
        }

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));
            String header = fileReader.readLine();

            final String[] columns = header.split(DELIMITER,-1);

            int alleleIndex = -1;
            int peptideIndex = -1;

            // all other fields are optional
            Integer upFlankIndex = null;
            Integer downFlankIndex = null;
            Integer sourceIndex = null;
            Integer tpmIndex = null;

            int otherDataCount = 0;
            List<Integer> otherColumnIndices = Lists.newArrayList();

            for(int i = 0; i < columns.length; ++i)
            {
                String column = columns[i];

                if(column.equals(FLD_ALLELE))
                    alleleIndex = i;
                else if(column.equals(FLD_PEPTIDE))
                    peptideIndex = i;
                else if(column.equals(FLD_UP_FLANK))
                    upFlankIndex = i;
                else if(column.equals(FLD_DOWN_FLANK))
                    downFlankIndex = i;
                else if(column.equals(FLD_SOURCE))
                    sourceIndex = i;
                else if(column.equals(FLD_TPM))
                    tpmIndex = i;
                else
                {
                    otherColumnIndices.add(i);
                    otherColumns.put(column, otherDataCount++);
                }
            }

            int itemCount = 0;

            String line = "";

            while((line = fileReader.readLine()) != null)
            {
                final String[] values = line.split(DELIMITER, -1);

                String allele = cleanAllele(values[alleleIndex]);
                String peptide = values[peptideIndex];

                if(!restrictedLengths.isEmpty() && !restrictedLengths.contains(peptide.length()))
                    continue;

                if(peptide.contains(AMINO_ACID_21ST)) // or stop codon
                    continue;

                ++itemCount;

                String source = sourceIndex != null ? values[sourceIndex] : "";
                String upFlank = upFlankIndex != null ? values[upFlankIndex] : "";
                String downFlank = downFlankIndex != null ? values[downFlankIndex] : "";

                BindData bindData = new BindData(allele, peptide, source, upFlank, downFlank);

                if(tpmIndex != null)
                {
                    bindData.setTPM(Double.parseDouble(values[tpmIndex]));
                }

                for(int i = 0; i < values.length; ++i)
                {
                    if(otherColumnIndices.contains(i))
                        bindData.getOtherData().add(values[i]);
                }

                Map<Integer,List<BindData>> pepLenBindDataMap = allelePeptideMap.get(allele);

                if(pepLenBindDataMap == null)
                {
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
