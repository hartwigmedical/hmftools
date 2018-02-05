package com.hartwig.hmftools.patientreporter.report.data;

import static java.util.stream.Collectors.collectingAndThen;
import static java.util.stream.Collectors.groupingBy;
import static java.util.stream.Collectors.joining;
import static java.util.stream.Collectors.mapping;
import static java.util.stream.Collectors.toList;
import static java.util.stream.Collectors.toSet;

import static net.sf.dynamicreports.report.builder.DynamicReports.field;

import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.apiclients.civic.data.CivicEvidenceItem;
import com.hartwig.hmftools.apiclients.civic.data.CivicVariant;
import com.hartwig.hmftools.common.gene.GeneCopyNumber;
import com.hartwig.hmftools.patientreporter.civic.AdditionalCivicMatches;
import com.hartwig.hmftools.patientreporter.variants.VariantReport;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import net.sf.dynamicreports.report.builder.FieldBuilder;

@Value.Immutable
@Value.Style(allParameters = true)
public abstract class Alteration {
    public static final FieldBuilder<?> ALTERATION = field("alteration", String.class);

    public abstract String getGene();

    public abstract String getPredictedEffect();

    public abstract List<AlterationEvidence> getEvidence();

    public abstract List<AlterationMatch> getMatches();

    @SuppressWarnings("unused")
    public String getAlteration() {
        // KODU: This function is used via reflection inside FieldBuilder -> pls do not remove!
        return getGene() + "\n" + getPredictedEffect();
    }

    @NotNull
    public static Alteration from(@NotNull final VariantReport variantReport, @NotNull final List<CivicVariant> civicVariants,
            @NotNull final Set<String> relevantDoids) {
        final List<AlterationEvidence> exactMatchEvidence = Lists.newArrayList();
        final List<AlterationMatch> matchingVariants = Lists.newArrayList();
        civicVariants.forEach(civicVariant -> {
            if (civicVariant.containsHgvsExpression(variantReport.hgvsProtein()) || civicVariant.coordinates()
                    .equals(variantReport.variant())) {
                exactMatchEvidence.addAll(alterationEvidence(civicVariant.evidenceItems(), civicVariant.summaryUrl(), relevantDoids));
                matchingVariants.add(AlterationMatch.of("exact", civicVariant));
            } else if (AdditionalCivicMatches.lossOfFunctionVariants(variantReport).contains(civicVariant.id())) {
                exactMatchEvidence.addAll(alterationEvidence(civicVariant.evidenceItems(), civicVariant.summaryUrl(), relevantDoids));
                matchingVariants.add(AlterationMatch.of("manual", civicVariant));
            } else if (civicVariant.coordinates().containsVariant(variantReport.variant())) {
                matchingVariants.add(AlterationMatch.of("approx.", civicVariant));
            }
        });
        return ImmutableAlteration.of(variantReport.gene(), variantReport.hgvsProtein(), exactMatchEvidence, matchingVariants);
    }

    @NotNull
    public static Alteration from(@NotNull final GeneCopyNumber copyNumberReport, @NotNull final List<CivicVariant> civicVariants,
            @NotNull final Set<String> relevantDoids) {
        final List<AlterationMatch> matchingVariants = Lists.newArrayList();
        final List<AlterationEvidence> exactMatchEvidence = Lists.newArrayList();
        civicVariants.forEach(civicVariant -> {
            if (AdditionalCivicMatches.copyNumberVariants(copyNumberReport).contains(civicVariant.id())) {
                exactMatchEvidence.addAll(alterationEvidence(civicVariant.evidenceItems(), civicVariant.summaryUrl(), relevantDoids));
                matchingVariants.add(AlterationMatch.of("manual", civicVariant));
            }
        });
        return ImmutableAlteration.of(copyNumberReport.gene(), copyNumberReport.alteration().name(), exactMatchEvidence, matchingVariants);
    }

    @NotNull
    private static List<AlterationEvidence> alterationEvidence(@NotNull final List<CivicEvidenceItem> evidenceItems,
            @NotNull final String civicUrl, @NotNull final Set<String> relevantDoids) {
        final Map<String, String> associationsPerSignificance = drugAssociationsPerSignificance(evidenceItems, relevantDoids);
        return associationsPerSignificance.keySet().stream().map(significance -> {
            final String associations = associationsPerSignificance.get(significance);
            return ImmutableAlterationEvidence.of(significance, associations, "CIViC", civicUrl);
        }).collect(toList());
    }

    @NotNull
    private static List<CivicEvidenceItem> onLabelEvidence(@NotNull final List<CivicEvidenceItem> evidenceItems,
            @NotNull Set<String> relevantDoids) {
        return evidenceItems.stream()
                .filter(item -> item.level() < 'D' && isUsable(item) && relevantDoids.contains(item.disease().doidString()))
                .sorted(Comparator.comparing(CivicEvidenceItem::level))
                .collect(toList());
    }

    @NotNull
    private static List<CivicEvidenceItem> offLabelEvidence(@NotNull final List<CivicEvidenceItem> evidenceItems,
            @NotNull Set<String> relevantDoids) {
        return evidenceItems.stream()
                .filter(item -> item.level() < 'C' && isUsable(item) && !relevantDoids.contains(item.disease().doidString()))
                .sorted(Comparator.comparing(CivicEvidenceItem::level))
                .collect(toList());
    }

    private static boolean isSupports(@Nullable final String evidenceDirection) {
        return evidenceDirection != null && evidenceDirection.toLowerCase().equals("supports");
    }

    private static boolean hasSignificance(@Nullable final String evidenceSignificance) {
        return evidenceSignificance != null && !evidenceSignificance.isEmpty() && !evidenceSignificance.toLowerCase().equals("n/a");
    }

    private static boolean isUsable(@NotNull final CivicEvidenceItem item) {
        return item.status().toLowerCase().equals("accepted") && isSupports(item.direction()) && hasSignificance(item.significance())
                && !item.drugs().isEmpty();
    }

    //MIVO: evidence item levels grouped by significance (sensitivity/resistance) and by drug.
    @NotNull
    private static Map<String, Map<String, String>> groupOnLabelEvidenceLevels(@NotNull final List<CivicEvidenceItem> evidenceItems) {
        return evidenceItems.stream()
                .flatMap(evidenceItem -> evidenceItem.drugNames().stream().map(drug -> ImmutablePair.of(drug, evidenceItem)))
                .collect(groupingBy(pair -> pair.getRight().significance(), groupingBy(Pair::getLeft,
                        mapping(pair -> pair.getRight().level(), collectingAndThen(toList(), Alteration::flattenOnLabelEvidenceLevels)))));
    }

    //MIVO: highest evidence level by significance (sensitivity/resistance) and by drug.
    @NotNull
    private static Map<String, Map<String, String>> groupOffLabelEvidenceLevels(@NotNull final List<CivicEvidenceItem> onLabelEvidenceItems,
            @NotNull final List<CivicEvidenceItem> offLabelEvidenceItems) {
        final Set<String> onLabelDrugs =
                onLabelEvidenceItems.stream().flatMap(evidenceItem -> evidenceItem.drugNames().stream()).collect(toSet());
        return offLabelEvidenceItems.stream()
                .flatMap(evidenceItem -> evidenceItem.drugNames()
                        .stream()
                        .filter(drug -> !onLabelDrugs.contains(drug))
                        .map(drug -> ImmutablePair.of(drug, evidenceItem)))
                .collect(groupingBy(pair -> pair.getRight().significance(), groupingBy(Pair::getLeft,
                        mapping(pair -> pair.getRight().level(), collectingAndThen(toList(), Alteration::flattenOffLabelEvidenceLevels)))));
    }

    @NotNull
    private static Map<String, String> drugAssociationsPerSignificance(@NotNull final List<CivicEvidenceItem> evidenceItems,
            @NotNull Set<String> relevantDoids) {
        final Map<String, String> results = Maps.newHashMap();
        final List<CivicEvidenceItem> onLabelEvidence = onLabelEvidence(evidenceItems, relevantDoids);
        final List<CivicEvidenceItem> offLabelEvidence = offLabelEvidence(evidenceItems, relevantDoids);
        final Map<String, Map<String, String>> onLabelEvidenceLevels = groupOnLabelEvidenceLevels(onLabelEvidence);
        final Map<String, Map<String, String>> offLabelEvidenceLevels = groupOffLabelEvidenceLevels(onLabelEvidence, offLabelEvidence);
        onLabelEvidenceLevels.forEach((significance, drugLevelsMap) -> results.put(significance, flattenDrugLevelsMap(drugLevelsMap)));
        offLabelEvidenceLevels.forEach((significance, drugLevelsMap) -> results.merge(significance, flattenDrugLevelsMap(drugLevelsMap),
                (existingAssociation, newAssociation) -> existingAssociation + "\n" + newAssociation));
        return results;
    }

    @NotNull
    private static String flattenOnLabelEvidenceLevels(@NotNull final List<Character> levels) {
        return "(" + levels.stream().distinct().map(String::valueOf).collect(joining(",")) + ")";
    }

    @NotNull
    private static String flattenOffLabelEvidenceLevels(@NotNull final List<Character> levels) {
        //MIVO: assumes levels list is sorted and non-empty
        return "(off-label: " + levels.get(0) + ")";
    }

    @NotNull
    private static String flattenDrugLevelsMap(@NotNull final Map<String, String> drugToLevels) {
        return drugToLevels.entrySet().stream().map(entry -> entry.getKey() + " " + entry.getValue()).collect(joining("\n"));
    }
}
