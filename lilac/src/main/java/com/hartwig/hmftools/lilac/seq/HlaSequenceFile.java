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
    public final List<HlaSequence> readFile(final String filename)
    {
        if(filename == null)
            return Lists.newArrayList();

        try
        {
            final List<String> fileData = Files.readAllLines(new File(filename).toPath());

            //fileData.remove(0);
            final Map<String, HlaSequence> entries = Maps.newHashMap();

            for(final String line : fileData) {
                String lineData = line.trim();

                if(lineData.startsWith("*", 1))
                    continue;

                String[] split = lineData.split(" ");
                String alleleStr = split[0].trim();
                int alleleIndex = lineData.indexOf(alleleStr);
                String remainder = lineData.substring(alleleIndex + alleleStr.length()).trim().replace(" ", "");

                HlaSequence sequence = entries.get(alleleStr);

                if(sequence == null)
                    entries.put(alleleStr, new HlaSequence(HlaAllele.fromString(alleleStr), remainder));
                else
                    sequence.appendSequence(remainder);
            }

            return entries.values().stream().collect(Collectors.toList());
        }
        catch (IOException e)
        {
            LL_LOGGER.error("failed to read HLF nucleotide file({}): {}", filename, e.toString());
            return Lists.newArrayList();
        }
    }

    @NotNull
    public final List<HlaSequence> reduceToFourDigit(final List<HlaSequence> sequences)
    {
        return null; // reduce($receiver, (Function1<? super HlaAllele, HlaAllele>) ((Function1) reduceToFourDigit .1.INSTANCE));
    }

    @NotNull
    public final List<HlaSequence> reduceToSixDigit(final List<HlaSequence> sequences)
    {
        return null; // this.reduce($receiver, (Function1<? super HlaAllele, HlaAllele>) ((Function1) reduceToSixDigit .1.INSTANCE));
    }

    /*
    private final List<HlaSequence> reduce(List<HlaSequence> $receiver, Function1<? super HlaAllele, HlaAllele> transform)
    {
        LinkedHashMap resultMap = new LinkedHashMap();
        for(HlaSequence sequence : $receiver)
        {
            HlaAllele reducedAllele = (HlaAllele) transform.invoke((Object) sequence.getAllele());
            if(resultMap.containsKey(reducedAllele))
            {
                continue;
            }
            ((Map) resultMap).put(reducedAllele, new HlaSequence(reducedAllele, sequence.getRawSequence()));
        }
        Collection collection = resultMap.values();

        return CollectionsKt.toList((Iterable) collection);
    }

     */

}
