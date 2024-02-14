package com.hartwig.hmftools.peach.data_loader;

import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.peach.PeachUtils.PCH_LOGGER;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.peach.PeachUtils;
import com.hartwig.hmftools.peach.effect.HaplotypeFunctionality;
import com.hartwig.hmftools.peach.effect.ImmutableHaplotypeFunctionality;

import org.jetbrains.annotations.NotNull;

public class HaplotypeFunctionalityLoader
{
    @NotNull
    public static List<HaplotypeFunctionality> loadFunctionalities(@NotNull String filename)
    {
        try
        {
            List<String> lines = Files.readAllLines(Paths.get(filename));
            List<HaplotypeFunctionality> functionalities = fromLines(lines);

            PCH_LOGGER.info("loaded {} lines of haplotype functionality from file ({})", functionalities.size(), filename);
            return functionalities;
        }
        catch(Exception e)
        {
            throw new RuntimeException(String.format("failed to load drugs TSV: %s", filename), e);
        }
    }

    @NotNull
    private static List<HaplotypeFunctionality> fromLines(@NotNull List<String> lines)
    {
        String header = lines.get(0);
        Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(header, PeachUtils.TSV_DELIMITER);

        List<HaplotypeFunctionality> functionalities = new ArrayList<>();
        for(String line : lines.subList(1, lines.size()))
        {
            String[] splitLine = line.split(PeachUtils.TSV_DELIMITER);
            HaplotypeFunctionality functionality = ImmutableHaplotypeFunctionality.builder()
                    .geneName(splitLine[fieldsIndexMap.get("Gene")])
                    .haplotypeName(splitLine[fieldsIndexMap.get("Haplotype")])
                    .functionality(splitLine[fieldsIndexMap.get("Functionality")])
                    .build();
            functionalities.add(functionality);
        }
        return functionalities;
    }
}
