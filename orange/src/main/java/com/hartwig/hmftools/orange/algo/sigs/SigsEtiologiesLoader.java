package com.hartwig.hmftools.orange.algo.sigs;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;

import org.jetbrains.annotations.NotNull;

public final class SigsEtiologiesLoader
{
    @NotNull
    public static Map<String, String> read(@NotNull String tsv) throws IOException
    {
        return fromLines(Files.readAllLines(new File(tsv).toPath()));
    }

    @NotNull
    private static Map<String, String> fromLines(@NotNull List<String> lines)
    {
        Map<String, String> configRules = Maps.newHashMap();

        Map<String, Integer> fields = createFieldsIndexMap(lines.get(0), TSV_DELIM);

        for(String line : lines.subList(1, lines.size()))
        {
            String[] values = line.split(TSV_DELIM, -1);

            configRules.put(values[fields.get("signature")], values[fields.get("etiology")]);
        }

        return configRules;
    }
}
