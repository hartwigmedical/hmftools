package com.hartwig.hmftools.lilac.seq;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.hla.HlaGene;

import org.jetbrains.annotations.NotNull;

public final class HlaSequenceFile
{
    public static final char SEQUENCE_DELIM = '|';

    private HlaSequenceFile() {}

    public static HlaSequenceLoci createFromReference(final HlaAllele allele, final String sequenceStr, final boolean isProteinFile)
    {
        List<String> sequences = isProteinFile ?
                Lists.newArrayListWithExpectedSize(367) : Lists.newArrayListWithExpectedSize(1098);

        int index = 0;
        StringBuilder sequence = new StringBuilder();
        boolean inMulti = false; // to handle inserts, ie more than 1 char, indicated by splitting the sequence by '|' chars
        while(index < sequenceStr.length())
        {
            char nextChar = sequenceStr.charAt(index);
            boolean isMulti = nextChar == SEQUENCE_DELIM;

            if(inMulti || isMulti)
            {
                if(inMulti && isMulti)
                {
                    inMulti = false;
                    sequences.add(sequence.toString());
                    sequence.setLength(0);
                }
                else if(isMulti)
                {
                    // start of new multi-char sequence
                    inMulti = true;
                }
                else
                {
                    sequence.append(nextChar);
                }
            }
            else
            {
                sequences.add(String.valueOf(nextChar));
            }

            ++index;
        }

        return new HlaSequenceLoci(allele, sequences);
    }

    public static List<HlaSequence> readDefintionFile(final File filename, final HlaGene gene)
    {
        if(filename == null)
            return Lists.newArrayList();

        try
        {
            String linePrefix = gene.shortName() + "*";
            final List<String> fileData = Files.readAllLines(filename.toPath());

            final List<String> orderedAlleles = Lists.newArrayList();
            final Map<String, HlaSequence> entries = Maps.newHashMap();

            for(final String line : fileData)
            {
                String lineData = line.trim();

                if(!lineData.startsWith(linePrefix))
                    continue;

                String[] split = lineData.split(" ");
                String alleleStr = split[0].trim();
                String remainder = lineData.substring(alleleStr.length()).trim().replace(" ", "");

                HlaSequence sequence = entries.get(alleleStr);

                if(sequence == null)
                {
                    orderedAlleles.add(alleleStr);
                    entries.put(alleleStr, new HlaSequence(HlaAllele.fromString(alleleStr), remainder));
                }
                else
                {
                    sequence.appendSequence(remainder);
                }
            }

            return orderedAlleles.stream().map(entries::get).collect(Collectors.toList());

            // return entries.values().stream().collect(Collectors.toList());
        }
        catch(IOException e)
        {
            LL_LOGGER.error("failed to read HLF nucleotide file({}): {}", filename.toString(), e.toString());
            return Lists.newArrayList();
        }
    }

    @NotNull
    public static List<HlaSequence> reduceToFourDigit(final Iterable<HlaSequence> sequences)
    {
        return reduce(sequences, true);
    }

    @NotNull
    public static List<HlaSequence> reduceToSixDigit(final Iterable<HlaSequence> sequences)
    {
        return reduce(sequences, false);
    }

    public static HlaAllele asSixDigit(final HlaAllele allele)
    {
        return new HlaAllele(allele.Gene, allele.AlleleGroup, allele.Protein, allele.Synonymous, "", null, null);
    }

    private static List<HlaSequence> reduce(final Iterable<HlaSequence> sequences, boolean toFourDigit)
    {
        Map<String, HlaSequence> reducedMap = Maps.newHashMap();
        List<String> orderAlleles = Lists.newArrayList();

        for(HlaSequence sequence : sequences)
        {
            HlaAllele reducedAllele = toFourDigit ? sequence.Allele.asFourDigit() : asSixDigit(sequence.Allele);

            if(reducedMap.containsKey(reducedAllele.toString()))
                continue;

            orderAlleles.add(reducedAllele.toString());
            reducedMap.put(reducedAllele.toString(), new HlaSequence(reducedAllele, sequence.getRawSequence()));
        }

        return orderAlleles.stream().map(reducedMap::get).collect(Collectors.toList());
    }
}
