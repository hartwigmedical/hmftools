package com.hartwig.hmftools.patientreporter.cfreport.data;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.actionability.WithEvent;
import com.hartwig.hmftools.common.actionability.variant.VariantEvidenceAnalyzer;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel;
import com.hartwig.hmftools.protect.variants.ReportableVariant;
import com.hartwig.hmftools.protect.variants.ReportableVariantSource;

import org.jetbrains.annotations.NotNull;

public final class EventFilter {

    private EventFilter() {
    }

    @NotNull
    public static <X extends WithEvent> List<X> removeEvidenceOnFilteredGermlineVariants(@NotNull List<X> evidences,
            @NotNull List<ReportableVariant> variants, @NotNull LimsGermlineReportingLevel germlineReportingLevel) {
        Set<String> germlineEvents = Sets.newHashSet();
        for (ReportableVariant variant : variants) {
            if (variant.source() == ReportableVariantSource.GERMLINE) {
                germlineEvents.add(VariantEvidenceAnalyzer.eventString(variant));
            }
        }

        List<X> filtered = Lists.newArrayList();
        for (X evidence : evidences) {
            if (germlineReportingLevel != LimsGermlineReportingLevel.NO_REPORTING || !germlineEvents.contains(evidence.event())) {
                filtered.add(evidence);
            }
        }

        return filtered;
    }
}
