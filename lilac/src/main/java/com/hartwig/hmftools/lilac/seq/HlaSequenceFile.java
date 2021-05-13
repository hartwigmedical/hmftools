package com.hartwig.hmftools.lilac.seq;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.lilac.hla.HlaAllele;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.jetbrains.annotations.NotNull;

public class HlaSequenceFile
{
    public static List<HlaSequence> readFile(final String filename)
    {
        if(filename == null)
            return Lists.newArrayList();

        try
        {
            final List<String> fileData = Files.readAllLines(new File(filename).toPath());

            final List<String> orderedAlleles = Lists.newArrayList();
            final Map<String,HlaSequence> entries = Maps.newHashMap();

            for(final String line : fileData)
            {
                String lineData = line.trim();

                if(!lineData.startsWith("*", 1))
                    continue;

                String[] split = lineData.split(" ");
                String alleleStr = split[0].trim();
                int alleleIndex = lineData.indexOf(alleleStr);
                String remainder = lineData.substring(alleleIndex + alleleStr.length()).trim().replace(" ", "");

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

            return orderedAlleles.stream().map(x -> entries.get(x)).collect(Collectors.toList());

            // return entries.values().stream().collect(Collectors.toList());
        }
        catch (IOException e)
        {
            LL_LOGGER.error("failed to read HLF nucleotide file({}): {}", filename, e.toString());
            return Lists.newArrayList();
        }
    }

    @NotNull
    public static List<HlaSequence> reduceToFourDigit(final List<HlaSequence> sequences)
    {
        return reduce(sequences, true);
    }

    @NotNull
    public static List<HlaSequence> reduceToSixDigit(final List<HlaSequence> sequences)
    {
        return reduce(sequences, false);
    }

    public static HlaAllele asSixDigit(final HlaAllele allele)
    {
        return new HlaAllele(allele.Gene, allele.AlleleGroup, allele.Protein, allele.Synonymous, "", null, null);
    }

    private static List<HlaSequence> reduce(final List<HlaSequence> sequences, boolean toFourDigit)
    {
        Map<String,HlaSequence> reducedMap = Maps.newHashMap();
        List<String> orderAlleles = Lists.newArrayList();

        for(HlaSequence sequence : sequences)
        {
            HlaAllele reducedAllele = toFourDigit ? sequence.Allele.asFourDigit() : asSixDigit(sequence.Allele);

            if(reducedMap.containsKey(reducedAllele.toString()))
                continue;

            orderAlleles.add(reducedAllele.toString());
            reducedMap.put(reducedAllele.toString(), new HlaSequence(reducedAllele, sequence.getRawSequence()));
        }

        return orderAlleles.stream().map(x -> reducedMap.get(x)).collect(Collectors.toList());
    }

}
