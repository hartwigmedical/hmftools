package com.hartwig.hmftools.patientreporter.cfreport.data;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.actionability.variant.VariantEvidenceAnalyzer;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel;
import com.hartwig.hmftools.protect.variants.ReportableVariant;
import com.hartwig.hmftools.protect.variants.ReportableVariantSource;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class EvidenceItems {

    private EvidenceItems() {
    }

    @NotNull
    public static List<EvidenceItem> removeEvidenceOnFilteredGermlineVariants(@NotNull List<EvidenceItem> evidenceItems,
            @NotNull List<ReportableVariant> variants, @NotNull LimsGermlineReportingLevel germlineReportingLevel) {
        Set<String> germlineEvents = Sets.newHashSet();
        for (ReportableVariant variant : variants) {
            if (variant.source() == ReportableVariantSource.GERMLINE) {
                germlineEvents.add(VariantEvidenceAnalyzer.eventString(variant));
            }
        }

        List<EvidenceItem> filtered = Lists.newArrayList();
        for (EvidenceItem item : evidenceItems) {
            if (germlineReportingLevel != LimsGermlineReportingLevel.NO_REPORTING || !germlineEvents.contains(item.event())) {
                filtered.add(item);
            }
        }

        return filtered;
    }

    @NotNull
    public static List<EvidenceItem> sort(@NotNull List<EvidenceItem> evidenceItems) {
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
    public static String sourceUrl(@NotNull EvidenceItem item) {
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

    public static int uniqueEventCount(@NotNull List<EvidenceItem> evidenceItems) {
        Set<String> events = Sets.newHashSet();
        for (EvidenceItem evidence : evidenceItems) {
            events.add(evidence.event());
        }
        return events.size();
    }

    public static int uniqueTherapyCount(@NotNull List<EvidenceItem> evidenceItems) {
        Set<String> drugs = Sets.newHashSet();
        for (EvidenceItem evidence : evidenceItems) {
            drugs.add(evidence.drug());
        }
        return drugs.size();
    }
}
