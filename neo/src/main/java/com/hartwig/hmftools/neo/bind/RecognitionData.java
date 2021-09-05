package com.hartwig.hmftools.neo.bind;

import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.DELIMITER;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindCommon.DELIM;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_ALLELE;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_IMMUNOGENIC;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_PEPTIDE;
import static com.hartwig.hmftools.neo.bind.BindCommon.cleanAllele;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;

public class RecognitionData
{
    public final String Allele;
    public final String Peptide;
    public final boolean Immunogenic;

    // optional fields
    private final List<String> mOtherData;

    public RecognitionData(final String allele, final String peptide, final boolean immunogenic)
    {
        Allele = allele;
        Peptide = peptide;
        Immunogenic = immunogenic;

        mOtherData = Lists.newArrayList();
    }

    public int peptideLength() { return Peptide.length(); }

    public final List<String> getOtherData() { return mOtherData; }

    public String toString()
    {
        return String.format("allele(%s) peptide(%s) immunogenic(%s)", Allele, Peptide, Immunogenic);
    }

    public static boolean loadRecognitionData(final String filename, final List<RecognitionData> recognitionData)
    {
        // final Map<String,Integer> otherColumns

        if(filename == null || !Files.exists(Paths.get(filename)))
        {
            NE_LOGGER.error("recognition data file({}) not found", filename);
            return false;
        }

        try
        {
            List<String> lines = Files.readAllLines(Paths.get(filename));
            String header = lines.get(0);
            lines.remove(0);

            Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, DELIM);

            int alleleIndex = fieldsIndexMap.get(FLD_ALLELE);
            int peptideIndex = fieldsIndexMap.get(FLD_PEPTIDE);
            int immunoIndex = fieldsIndexMap.get(FLD_IMMUNOGENIC);

            for(String line :lines)
            {
                final String[] values = line.split(DELIMITER, -1);

                String allele = cleanAllele(values[alleleIndex]);
                String peptide = values[peptideIndex];
                boolean immunogenic = values[immunoIndex].equals("1") || Boolean.parseBoolean(values[immunoIndex]);

                recognitionData.add(new RecognitionData(allele, peptide, immunogenic));
            }

            NE_LOGGER.info("loaded {} peptides with immunogenic status from file({})", lines.size(), filename);
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to read immunogenic data file: {}", e.toString());
            return false;
        }

        return true;
    }
}
