package com.hartwig.hmftools.peach.data_loader;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.peach.PeachUtils.PCH_LOGGER;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.peach.effect.DrugInfo;
import com.hartwig.hmftools.peach.effect.ImmutableDrugInfo;

import org.jetbrains.annotations.NotNull;

public class DrugInfoLoader
{
    @NotNull
    public static List<DrugInfo> loadDrugInfos(@NotNull String filename)
    {
        try
        {
            List<String> lines = Files.readAllLines(Paths.get(filename));
            List<DrugInfo> drugInfos = fromLines(lines);

            PCH_LOGGER.info("loaded {} lines of drug info from file ({})", drugInfos.size(), filename);
            return drugInfos;
        }
        catch(Exception e)
        {
            throw new RuntimeException(String.format("failed to load drugs TSV: %s", filename), e);
        }
    }

    @NotNull
    private static List<DrugInfo> fromLines(@NotNull List<String> lines)
    {
        String header = lines.get(0);
        Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);

        List<DrugInfo> drugInfos = new ArrayList<>();
        for(String line : lines.subList(1, lines.size()))
        {
            String[] splitLine = line.split(TSV_DELIM);
            DrugInfo drugInfo = ImmutableDrugInfo.builder()
                    .drugName(splitLine[fieldsIndexMap.get("drugName")])
                    .geneName(splitLine[fieldsIndexMap.get("gene")])
                    .prescriptionInfoUrl(splitLine[fieldsIndexMap.get("urlPrescriptionInfo")])
                    .build();
            drugInfos.add(drugInfo);
        }
        return drugInfos;
    }
}
