package com.hartwig.hmftools.peach.output;

import com.hartwig.hmftools.peach.HaplotypeAnalysis;
import com.hartwig.hmftools.peach.effect.DrugInfoStore;
import com.hartwig.hmftools.peach.effect.HaplotypeFunctionStore;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import static com.hartwig.hmftools.peach.PeachUtils.GERMLINE_TOTAL_COPY_NUMBER;
import static com.hartwig.hmftools.peach.PeachUtils.TSV_DELIMITER;

public class BestHaplotypeCombinationsFile
{
    private static final String UNKNOWN_ALLELE_STRING = "UNRESOLVED";
    private static final String DRUG_SEPARATOR = ";";

    public static void write(@NotNull String filePath, @NotNull Map<String, HaplotypeAnalysis> geneToHaplotypeAnalysis,
            @Nullable DrugInfoStore drugInfoStore, @Nullable HaplotypeFunctionStore haplotypeFunctionStore) throws IOException
    {
        Files.write(new File(filePath).toPath(), toLines(geneToHaplotypeAnalysis, drugInfoStore, haplotypeFunctionStore));
    }

    @NotNull
    public static List<String> toLines(@NotNull Map<String, HaplotypeAnalysis> geneToHaplotypeAnalysis,
            @Nullable DrugInfoStore drugInfoStore, @Nullable HaplotypeFunctionStore haplotypeFunctionStore)
    {
        List<String> lines = new ArrayList<>();
        lines.add(header());
        geneToHaplotypeAnalysis.entrySet()
                .stream()
                .sorted(Map.Entry.comparingByKey())
                .map(e -> toLines(e.getKey(), e.getValue(), drugInfoStore, haplotypeFunctionStore))
                .flatMap(Collection::stream)
                .forEach(lines::add);
        return lines;
    }

    @NotNull
    private static String header()
    {
        return new StringJoiner(TSV_DELIMITER).add("gene")
                .add("allele")
                .add("count")
                .add("function")
                .add("linkedDrugs")
                .add("prescriptionUrls")
                .toString();
    }

    @NotNull
    private static List<String> toLines(@NotNull String gene, @NotNull HaplotypeAnalysis analysis, @Nullable DrugInfoStore drugInfoStore,
            @Nullable HaplotypeFunctionStore haplotypeFunctionStore)
    {
        if(analysis.hasBestHaplotypeCombination())
        {
            return analysis.getBestHaplotypeCombination()
                    .getHaplotypeNameToCount()
                    .entrySet()
                    .stream()
                    .sorted(Map.Entry.comparingByKey())
                    .map(e -> toLine(gene, e.getKey(), e.getValue(), drugInfoStore, haplotypeFunctionStore))
                    .collect(Collectors.toList());
        }
        else
        {
            return List.of(toLine(gene, UNKNOWN_ALLELE_STRING, GERMLINE_TOTAL_COPY_NUMBER, drugInfoStore, haplotypeFunctionStore));
        }
    }

    @NotNull
    private static String toLine(@NotNull String gene, @NotNull String haplotypeName, int count, @Nullable DrugInfoStore drugInfoStore,
            @Nullable HaplotypeFunctionStore haplotypeFunctionStore)
    {
        String printableHaplotypeFunction = getPrintableHaplotypeFunction(gene, haplotypeName, haplotypeFunctionStore);
        String printableDrugNamesString = getPrintableDrugNamesString(gene, drugInfoStore);
        String printablePresciptionUrlsString = getPrintablePresciptionUrlsString(gene, drugInfoStore);
        return new StringJoiner(TSV_DELIMITER).add(gene)
                .add(haplotypeName)
                .add(Integer.toString(count))
                .add(printableHaplotypeFunction)
                .add(printableDrugNamesString)
                .add(printablePresciptionUrlsString)
                .toString();
    }

    @NotNull
    private static String getPrintableHaplotypeFunction(@NotNull String gene, @NotNull String haplotypeName,
            @Nullable HaplotypeFunctionStore haplotypeFunctionStore)
    {
        if(haplotypeFunctionStore == null)
        {
            return Strings.EMPTY;
        }
        String functionality = haplotypeFunctionStore.getFunction(gene, haplotypeName);
        return Objects.requireNonNullElse(functionality, Strings.EMPTY);
    }

    @NotNull
    private static String getPrintableDrugNamesString(@NotNull String gene, @Nullable DrugInfoStore drugInfoStore)
    {
        if(drugInfoStore == null)
        {
            return Strings.EMPTY;
        }
        List<String> sortedLinkedDrugNames = getSortedLinkedDrugNames(gene, drugInfoStore);
        if(sortedLinkedDrugNames.isEmpty())
        {
            return Strings.EMPTY;
        }
        else
        {
            return String.join(DRUG_SEPARATOR, sortedLinkedDrugNames);
        }
    }

    @NotNull
    private static String getPrintablePresciptionUrlsString(@NotNull String gene, @Nullable DrugInfoStore drugInfoStore)
    {
        if(drugInfoStore == null)
        {
            return Strings.EMPTY;
        }
        List<String> sortedLinkedDrugNames = getSortedLinkedDrugNames(gene, drugInfoStore);
        if(sortedLinkedDrugNames.isEmpty())
        {
            return Strings.EMPTY;
        }
        else
        {
            return sortedLinkedDrugNames.stream()
                    .map(d -> drugInfoStore.getPrescriptionInfoUrl(gene, d))
                    .collect(Collectors.joining(DRUG_SEPARATOR));
        }
    }

    @NotNull
    private static List<String> getSortedLinkedDrugNames(@NotNull String gene, @NotNull DrugInfoStore drugInfoStore)
    {
        return drugInfoStore.getRelevantDrugNames(gene).stream().sorted().collect(Collectors.toList());
    }
}
