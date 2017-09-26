package com.hartwig.hmftools.patientreporter.report.data;

import static java.util.stream.Collectors.groupingBy;
import static java.util.stream.Collectors.mapping;
import static java.util.stream.Collectors.toList;

import static net.sf.dynamicreports.report.builder.DynamicReports.field;

import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.apiclients.civic.data.CivicEvidenceItem;
import com.hartwig.hmftools.apiclients.civic.data.CivicVariant;
import com.hartwig.hmftools.patientreporter.variants.VariantReport;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import net.sf.dynamicreports.report.builder.FieldBuilder;

@Value.Immutable
@Value.Style(allParameters = true)
public abstract class Alteration {
    private static final Logger LOGGER = LogManager.getLogger(Alteration.class);

    public static final FieldBuilder<?> ALTERATION = field("alteration", String.class);

    public abstract String getGene();

    public abstract String getPredictedEffect();

    public abstract List<AlterationEvidence> getEvidence();

    public abstract List<AlterationMatch> getMatches();

    public String getAlteration() {
        return getGene() + "\n" + getPredictedEffect();
    }

    public static Alteration from(@NotNull final VariantReport variantReport, @NotNull final List<CivicVariant> civicVariants,
            @NotNull final Set<String> tumorSubtypesDoids) {
        final String gene = variantReport.gene();
        final String predictedEffect = variantReport.hgvsProtein();
        final List<AlterationEvidence> exactMatchEvidence = Lists.newArrayList();
        final List<AlterationMatch> matchingVariants = Lists.newArrayList();

        civicVariants.forEach(civicVariant -> {
            if (civicVariant.coordinates().equals(variantReport.variant())) {
                final List<CivicEvidenceItem> relevantEvidence = relevantEvidenceForTumor(civicVariant.evidenceItems(), tumorSubtypesDoids);
                final Map<String, Map<String, List<CivicEvidenceItem>>> groupedEvidenceItems = groupEvidenceItems(relevantEvidence);
                for (final String significance : groupedEvidenceItems.keySet()) {
                    final String drugsString = getDrugsWithEvidenceLevel(groupedEvidenceItems.get(significance));
                    exactMatchEvidence.add(ImmutableAlterationEvidence.of(significance, drugsString, "CIViC"));
                }
                matchingVariants.add(AlterationMatch.of("exact", civicVariant));
            } else {
                matchingVariants.add(AlterationMatch.of("approx.", civicVariant));
            }
        });
        return ImmutableAlteration.of(gene, predictedEffect, exactMatchEvidence, matchingVariants);
    }

    @NotNull
    private static String getDrugsWithEvidenceLevel(@NotNull final Map<String, List<CivicEvidenceItem>> drugToEvidenceItems) {
        final List<String> drugs = drugToEvidenceItems.entrySet()
                .stream()
                .map(entry -> entry.getKey() + "(" + Strings.join(
                        entry.getValue().stream().map(CivicEvidenceItem::level).distinct().collect(Collectors.toList()), ',') + ")")
                .collect(Collectors.toList());
        return Strings.join(drugs, '\n');
    }

    @NotNull
    private static List<CivicEvidenceItem> relevantEvidenceForTumor(@NotNull final List<CivicEvidenceItem> evidenceItems,
            @NotNull Set<String> tumorSubtypesDoids) {
        return evidenceItems.stream()
                .filter(item -> item.level() < 'C' && isUsable(item) && tumorSubtypesDoids.contains(item.disease().doidString()))
                .sorted(Comparator.comparing(CivicEvidenceItem::level))
                .collect(Collectors.toList());
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

    //MIVO: group evidence items by significance (sensitivity/resistance) and by drug.
    @NotNull
    private static Map<String, Map<String, List<CivicEvidenceItem>>> groupEvidenceItems(
            @NotNull final List<CivicEvidenceItem> evidenceItems) {
        return evidenceItems.stream()
                .flatMap(evidenceItem -> evidenceItem.drugs().stream().map(drug -> ImmutablePair.of(drug.name(), evidenceItem)))
                .collect(groupingBy(pair -> pair.getRight().significance(), groupingBy(Pair::getLeft, mapping(Pair::getRight, toList()))));
    }
}
