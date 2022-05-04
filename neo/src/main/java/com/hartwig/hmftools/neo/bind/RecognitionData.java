package com.hartwig.hmftools.neo.bind;

import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.DELIMITER;
import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindCommon.DELIM;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_ALLELE;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_IMMUNOGENIC;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_PEPTIDE;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_SOURCE;
import static com.hartwig.hmftools.neo.bind.BindCommon.cleanAllele;
import static com.hartwig.hmftools.neo.bind.BindConstants.DEFAULT_PEPTIDE_LENGTHS;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

public class RecognitionData
{
    public final String Allele;
    public final String Peptide;
    public final String Source;
    public final boolean Immunogenic;

    public RecognitionData(final String allele, final String peptide, final String source, boolean immunogenic)
    {
        Allele = allele;
        Peptide = peptide;
        Immunogenic = immunogenic;
        Source = source;
    }

    public int peptideLength() { return Peptide.length(); }

    public String toString()
    {
        return String.format("allele(%s) peptide(%s) source(%s) immunogenic(%s)",
                Allele, Peptide, Source, Immunogenic);
    }

    public static boolean loadRecognitionData(final String filename, final List<RecognitionData> recognitionData)
    {
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
            Integer immunoIndex = fieldsIndexMap.get(FLD_IMMUNOGENIC);
            Integer sourceIndex = fieldsIndexMap.get(FLD_SOURCE);

            for(String line :lines)
            {
                final String[] values = line.split(DELIMITER, -1);

                String allele = cleanAllele(values[alleleIndex]);
                String peptide = values[peptideIndex];

                if(!DEFAULT_PEPTIDE_LENGTHS.contains(peptide.length()))
                    continue;

                boolean immunogenic = immunoIndex != null ?
                        (values[immunoIndex].equals("1") || Boolean.parseBoolean(values[immunoIndex])) : true;

                String source = sourceIndex != null ? values[sourceIndex] : "";

                recognitionData.add(new RecognitionData(allele, peptide, source, immunogenic));
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
