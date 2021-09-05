package com.hartwig.hmftools.neo.bind;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindConstants.AMINO_ACIDS;
import static com.hartwig.hmftools.neo.bind.BindConstants.AMINO_ACID_COUNT;
import static com.hartwig.hmftools.neo.bind.BindConstants.INVALID_AMINO_ACID;
import static com.hartwig.hmftools.neo.bind.BindConstants.aminoAcidIndex;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Map;

import com.google.common.collect.Maps;

public class RecognitionCounts
{
    private final Map<String,Map<Integer,Map<Boolean,ImmuneCounts>>> mAlleleCounts; // counts data by peptide length>


    public RecognitionCounts()
    {
        mAlleleCounts = Maps.newHashMap();
    }

    public void process(final RecognitionData recogData)
    {
        Map<Integer,Map<Boolean,ImmuneCounts>> pepLenCounts = mAlleleCounts.get(recogData.peptideLength());

        if(pepLenCounts == null)
        {
            pepLenCounts = Maps.newHashMap();
            mAlleleCounts.put(recogData.Allele, pepLenCounts);
        }

        Map<Boolean,ImmuneCounts> immuneTypes = pepLenCounts.get(recogData.peptideLength());

        if(immuneTypes == null)
        {
            immuneTypes = Maps.newHashMap();
            pepLenCounts.put(recogData.peptideLength(), immuneTypes);
        }

        ImmuneCounts immuneCounts = immuneTypes.get(recogData.Immunogenic);

        if(immuneCounts == null)
        {
            immuneCounts = new ImmuneCounts(recogData.Immunogenic, recogData.peptideLength());
            immuneTypes.put(recogData.Immunogenic, immuneCounts);
        }

        for(int pos = 0; pos < recogData.Peptide.length(); ++pos)
        {
            char aminoAcid = recogData.Peptide.charAt(pos);
            int aaIndex = aminoAcidIndex(aminoAcid);

            if(aaIndex == INVALID_AMINO_ACID)
                continue;

            ++immuneCounts.Counts[aaIndex][pos];
        }
    }

    public void writeCounts(final String filename, int maxPeptideLength)
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("Allele,PeptideLength,Immunogenic,AminoAcid");

            for(int i = 0; i < maxPeptideLength; ++i)
            {
                writer.write(String.format(",P%d", i));
            }

            writer.newLine();

            for(Map.Entry<String,Map<Integer,Map<Boolean,ImmuneCounts>>> alleleEntry : mAlleleCounts.entrySet())
            {
                String allele = alleleEntry.getKey();

                for(Map<Boolean,ImmuneCounts> immuneCountsMap : alleleEntry.getValue().values())
                {
                    for(ImmuneCounts immuneCounts : immuneCountsMap.values())
                    {
                        for(int aa = 0; aa < AMINO_ACID_COUNT; ++aa)
                        {
                            char aminoAcid = AMINO_ACIDS.get(aa);

                            writer.write(String.format("%s,%d,%c", allele, immuneCounts.PeptideLength, aminoAcid));

                            for(int pos = 0; pos < maxPeptideLength; ++pos)
                            {
                                if(pos < immuneCounts.PeptideLength)
                                {
                                    writer.write(String.format(",%d", immuneCounts.Counts[aa][pos]));
                                }
                                else
                                {
                                    writer.write(",0");
                                }
                            }

                            writer.newLine();
                        }
                    }
                }
            }


        }
        catch (IOException e)
        {
            NE_LOGGER.error("failed to write immunogenic counts data: {}", e.toString());
        }
    }

    private class ImmuneCounts
    {
        public final boolean Immunogenic;
        public final int PeptideLength;
        public final int[][] Counts;

        public ImmuneCounts(final boolean immunogenic, final int peptideLength)
        {
            Immunogenic = immunogenic;
            PeptideLength = peptideLength;
            Counts = new int[AMINO_ACID_COUNT][PeptideLength];
        }
    }
}
