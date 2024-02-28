package com.hartwig.hmftools.peach.data_loader;

import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.peach.PeachUtils.PCH_LOGGER;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.peach.PeachUtils;
import com.hartwig.hmftools.peach.effect.HaplotypeFunction;
import com.hartwig.hmftools.peach.effect.ImmutableHaplotypeFunction;

import org.jetbrains.annotations.NotNull;

public class HaplotypeFunctionLoader
{
    @NotNull
    public static List<HaplotypeFunction> loadFunctions(@NotNull String filename)
    {
        try
        {
            List<String> lines = Files.readAllLines(Paths.get(filename));
            List<HaplotypeFunction> functions = fromLines(lines);

            PCH_LOGGER.info("loaded {} lines of haplotype functions from file ({})", functions.size(), filename);
            return functions;
        }
        catch(Exception e)
        {
            throw new RuntimeException(String.format("failed to load drugs TSV: %s", filename), e);
        }
    }

    @NotNull
    private static List<HaplotypeFunction> fromLines(@NotNull List<String> lines)
    {
        String header = lines.get(0);
        Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(header, PeachUtils.TSV_DELIMITER);

        List<HaplotypeFunction> functions = new ArrayList<>();
        for(String line : lines.subList(1, lines.size()))
        {
            String[] splitLine = line.split(PeachUtils.TSV_DELIMITER);
            HaplotypeFunction function = ImmutableHaplotypeFunction.builder()
                    .geneName(splitLine[fieldsIndexMap.get("Gene")])
                    .haplotypeName(splitLine[fieldsIndexMap.get("Haplotype")])
                    .function(splitLine[fieldsIndexMap.get("Function")])
                    .build();
            functions.add(function);
        }
        return functions;
    }
}
