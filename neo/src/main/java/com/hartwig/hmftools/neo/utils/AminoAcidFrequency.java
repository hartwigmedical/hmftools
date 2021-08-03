package com.hartwig.hmftools.neo.utils;

import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.DELIMITER;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.neo.bind.BindConstants;

public class AminoAcidFrequency
{
    private final Map<Character,Double> mAminoAcidFrequencies;

    public AminoAcidFrequency()
    {
        mAminoAcidFrequencies = Maps.newHashMap();
        loadFrequencies();
    }

    public Map<Character,Double> getAminoAcidFrequencies() { return mAminoAcidFrequencies; }

    public double getAminoAcidFrequency(final char aminoAcid)
    {
        Double percent = mAminoAcidFrequencies.get(aminoAcid);
        return percent != null ? percent : 1.0 / BindConstants.AMINO_ACIDS.size();
    }

    private void loadFrequencies()
    {
        final List<String> lines = new BufferedReader(new InputStreamReader(
                RefGenomeCoordinates.class.getResourceAsStream("/ref/amino_acid_frequencies.csv")))
                .lines().collect(Collectors.toList());

        final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), DELIMITER);
        lines.remove(0);

        int aaIndex = fieldsIndexMap.get("AminoAcid");
        int percentIndex = fieldsIndexMap.get("Percent");

        for(String line : lines)
        {
            final String[] items = line.split(DELIMITER, -1);

            char aminoAcid = items[aaIndex].charAt(0);
            double percent = Double.parseDouble(items[percentIndex]);
            mAminoAcidFrequencies.put(aminoAcid, percent);
        }
    }
}
