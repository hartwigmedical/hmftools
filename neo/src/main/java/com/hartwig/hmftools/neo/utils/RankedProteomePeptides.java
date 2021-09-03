package com.hartwig.hmftools.neo.utils;

import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.DELIMITER;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindCommon.DELIM;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_ALLELE;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_LIKE_RANK;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_PEPTIDE;
import static com.hartwig.hmftools.neo.bind.BindCommon.cleanAllele;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Map;

import com.google.common.collect.Maps;

public class RankedProteomePeptides
{
    private final Map<String,Map<Integer,Map<String,Double>>> mAllelePeptideRanks;
    private boolean mIsValid;

    public static final String PROTEOME_RANKS_FILE = "proteome_ranks_file";

    public RankedProteomePeptides(final String filename)
    {
        mAllelePeptideRanks = Maps.newHashMap();
        mIsValid = load(filename);
    }

    public boolean isValid() { return mIsValid; }

    public Map<String,Double> getPeptideRanks(final String allele, int peptideLength)
    {
        Map<Integer,Map<String,Double>> pepLenBindDataMap = mAllelePeptideRanks.get(allele);

        if(pepLenBindDataMap == null)
            return null;

        return pepLenBindDataMap.get(peptideLength);
    }

    private boolean load(final String filename)
    {
        if(filename == null || !Files.exists(Paths.get(filename)))
        {
            NE_LOGGER.error("proteome ranks file({}) not found", filename);
            return false;
        }

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));
            String header = fileReader.readLine();

            Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, DELIM);

            int alleleIndex = fieldsIndexMap.get(FLD_ALLELE);
            int peptideIndex = fieldsIndexMap.get(FLD_PEPTIDE);
            int rankIndex = fieldsIndexMap.get(FLD_LIKE_RANK);

            String line = "";

            while((line = fileReader.readLine()) != null)
            {
                final String[] values = line.split(DELIMITER, -1);

                String allele = cleanAllele(values[alleleIndex]);
                String peptide = values[peptideIndex];
                int peptideLength = peptide.length();
                double rank = Double.parseDouble(values[rankIndex]);

                Map<Integer,Map<String,Double>> pepLenBindDataMap = mAllelePeptideRanks.get(allele);

                if(pepLenBindDataMap == null)
                {
                    pepLenBindDataMap = Maps.newHashMap();
                    mAllelePeptideRanks.put(allele, pepLenBindDataMap);
                }

                Map<String,Double> peptideMap = pepLenBindDataMap.get(peptideLength);

                if(peptideMap == null)
                {
                    peptideMap = Maps.newHashMap();
                    pepLenBindDataMap.put(peptideLength, peptideMap);
                }

                peptideMap.put(peptide, rank);
            }

            NE_LOGGER.info("loaded {} alleles with proteome ranks from file({})", mAllelePeptideRanks.size(), filename);
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to read allele proteome rank file: {}", e.toString());
            return false;
        }

        return true;
    }
}
