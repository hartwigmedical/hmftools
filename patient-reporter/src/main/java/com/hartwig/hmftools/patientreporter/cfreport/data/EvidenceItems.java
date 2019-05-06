package com.hartwig.hmftools.patientreporter.cfreport.data;

import com.hartwig.hmftools.common.actionability.EvidenceItem;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.stream.Collectors;

public final class EvidenceItems {

    /**
     * Get all evidence items from the list that are not a trial source or should be included in report
     */
    @NotNull
    public static List<EvidenceItem> filter(@NotNull final List<EvidenceItem> evidenceItems) {
        return evidenceItems.stream()
                .filter(evidenceItem -> (!evidenceItem.source().isTrialSource() && evidenceItem.level().includeInReport()))
                .collect(Collectors.toList());
    }

    @NotNull
    public static List<EvidenceItem> sort(@NotNull final List<EvidenceItem> evidenceItems) {
        return evidenceItems.stream().sorted((item1, item2) -> {
            if (item1.level().equals(item2.level())) {
                if (item1.event().equals(item2.event())) {
                    return item1.drug().compareTo(item2.drug());
                } else {
                    return item1.event().compareTo(item2.event());
                }
            } else {
                return item1.level().readableString().compareTo(item2.level().readableString());
            }
        }).collect(Collectors.toList());
    }

    @NotNull
    public static String sourceUrl(@NotNull final EvidenceItem item) {
        String source = item.source().sourceName();
        String reference = item.reference();
        String gene = item.event();
        switch (source.toLowerCase()) {
            case "oncokb":
                String[] geneId = gene.split(" ");
                String referenceFormatting = reference.replace(" ", "%20");
                return "http://oncokb.org/#/gene/" + geneId[0] + "/alteration/" + referenceFormatting;
            case "cgi":
                return "https://www.cancergenomeinterpreter.org/biomarkers";
            case "civic":
                String[] variantId = reference.split(":");
                return "https://civic.genome.wustl.edu/links/variants/" + variantId[1];
            default:
                return Strings.EMPTY;
        }

    }

    public static int uniqueEventCount(@NotNull final List<EvidenceItem> evidenceItems) {
        return (int) evidenceItems.stream().filter(distinctByKey(EvidenceItem::event)).count();
    }

    public static int uniqueTherapyCount(@NotNull final List<EvidenceItem> evidenceItems) {
        return (int) evidenceItems.stream().filter(distinctByKey(EvidenceItem::drug)).count();
    }

    private static <T> Predicate<T> distinctByKey(@NotNull Function<? super T, Object> keyExtractor) {
        Map<Object, Boolean> map = new ConcurrentHashMap<>();
        return t -> map.putIfAbsent(keyExtractor.apply(t), Boolean.TRUE) == null;
    }
}
