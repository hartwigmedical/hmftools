package com.hartwig.hmftools.neo.utils;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.CSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.neo.bind.BindConstants;
import static com.hartwig.hmftools.neo.bind.BindConstants.AMINO_ACIDS;

import com.google.common.collect.Lists;

public class AminoAcidFrequency
{
    private final Map<Character,Double> mAminoAcidFrequencyMap;
    private final List<Double> mAminoAcidFrequencies; // indexed by amino acid from constants

    public AminoAcidFrequency()
    {
        mAminoAcidFrequencyMap = Maps.newHashMap();
        mAminoAcidFrequencies = Lists.newArrayList();
        loadFrequencies();
    }

    public Map<Character,Double> getAminoAcidFrequencies() { return mAminoAcidFrequencyMap; }

    public double getAminoAcidFrequency(final char aminoAcid)
    {
        Double percent = mAminoAcidFrequencyMap.get(aminoAcid);
        return percent != null ? percent : 1.0 / BindConstants.AMINO_ACIDS.size();
    }

    public double getAminoAcidFrequency(final int aaIndex)
    {
        return mAminoAcidFrequencies.get(aaIndex);
    }

    private void loadFrequencies()
    {
        final List<String> lines = new BufferedReader(new InputStreamReader(
                AminoAcidFrequency.class.getResourceAsStream("/ref/amino_acid_frequencies.csv")))
                .lines().collect(Collectors.toList());

        final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), CSV_DELIM);
        lines.remove(0);

        int aminoAcidIndex = fieldsIndexMap.get("AminoAcid");
        int percentIndex = fieldsIndexMap.get("Percent");

        for(String line : lines)
        {
            final String[] items = line.split(CSV_DELIM, -1);

            char aminoAcid = items[aminoAcidIndex].charAt(0);
            double percent = Double.parseDouble(items[percentIndex]);
            mAminoAcidFrequencyMap.put(aminoAcid, percent);
        }

        for(int aaIndex = 0 ; aaIndex < AMINO_ACIDS.size(); ++aaIndex)
        {
            char aminoAcid = AMINO_ACIDS.get(aaIndex);
            mAminoAcidFrequencies.add(mAminoAcidFrequencyMap.get(aminoAcid));
        }
    }
}
