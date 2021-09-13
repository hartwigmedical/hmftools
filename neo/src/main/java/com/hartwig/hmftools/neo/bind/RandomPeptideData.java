package com.hartwig.hmftools.neo.bind;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindCommon.DELIM;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_DOWN_FLANK;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_PEPTIDE;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_TPM;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_UP_FLANK;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;

public class RandomPeptideData
{
    public final String Peptide;
    public final String UpFlank;
    public final String DownFlank;
    public final double Expression;

    public RandomPeptideData(final String peptide, final String upFlank, final String downFlank, double expression)
    {
        Peptide = peptide;
        UpFlank = upFlank;
        DownFlank = downFlank;
        Expression = expression;
    }

    public static boolean loadRandomPeptides(final String filename, final Set<Integer> peptideLengths, final Map<Integer,List<RandomPeptideData>> randomPeptideMap)
    {
        if(!randomPeptideMap.isEmpty())
            return true;

        if(filename == null)
        {
            NE_LOGGER.error("missing random peptides file");
            return false;
        }

        try
        {
            List<String> lines = Files.readAllLines(new File(filename).toPath());

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), DELIM);
            lines.remove(0);

            int peptideIndex = fieldsIndexMap.get(FLD_PEPTIDE);
            Integer upFlankIndex = fieldsIndexMap.get(FLD_UP_FLANK);
            Integer downFlankIndex = fieldsIndexMap.get(FLD_DOWN_FLANK);
            Integer tpmIndex = fieldsIndexMap.get(FLD_TPM);

            if(upFlankIndex != null && downFlankIndex != null)
            {
                List<RandomPeptideData> peptideList = null;

                for(String line : lines)
                {
                    String[] values = line.split(DELIM, -1);
                    String peptide = values[peptideIndex];

                    int peptideLength = peptide.length();
                    peptideList = randomPeptideMap.get(peptideLength);

                    if(peptideList == null)
                    {
                        peptideList = Lists.newArrayList();
                        randomPeptideMap.put(peptideLength, peptideList);
                    }

                    double tpm = tpmIndex != null ? Double.parseDouble(values[tpmIndex]) : 0;

                    peptideList.add(new RandomPeptideData(
                            peptide, values[upFlankIndex], values[downFlankIndex], tpm));
                }
            }
            else
            {
                for(Integer peptideLength : peptideLengths)
                {
                    List<RandomPeptideData> peptideList = Lists.newArrayList();
                    randomPeptideMap.put(peptideLength, peptideList);

                    for(String line : lines)
                    {
                        String[] values = line.split(DELIM, -1);
                        String peptide = values[peptideIndex].substring(0, peptideLength);
                        peptideList.add(new RandomPeptideData(peptide, "", "", 0));
                    }
                }
            }

            NE_LOGGER.info("loaded {} random peptides", lines.size());
            return true;
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to load random peptide file: {}", e.toString());
            return false;
        }
    }

}
