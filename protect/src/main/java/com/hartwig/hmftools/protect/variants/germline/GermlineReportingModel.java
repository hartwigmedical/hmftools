package com.hartwig.hmftools.protect.variants.germline;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class GermlineReportingModel {

    private static final Logger LOGGER = LogManager.getLogger(GermlineReportingModel.class);

    @NotNull
    private final Set<String> reportableGenes = Sets.newHashSet();
    private final Set<String> notifiableGenes = Sets.newHashSet();
    private final Set<String> monoallelicGenes = Sets.newHashSet();
    private final Map<String, String> reportableSpecificVariants = Maps.newHashMap();

    public GermlineReportingModel(@NotNull final List<GermlineReporting> germlineReportingEntries) {
        for (GermlineReporting entry : germlineReportingEntries) {
            reportableGenes.add(entry.gene());
            if (entry.notifyClinicalGeneticist()) {
                notifiableGenes.add(entry.gene());
            }
            if (!entry.reportBiallelicOnly()) {
                monoallelicGenes.add(entry.gene());
            }
            if (entry.reportableSpecificVariant() != null) {
                reportableSpecificVariants.put(entry.gene(), entry.reportableSpecificVariant());
            }
        }
    }

    @NotNull
    public Map<String, String> reportableSpecificVariants() {
        return reportableSpecificVariants;
    }

    @NotNull
    public Set<String> monoallelicGenesReportable() {
        return monoallelicGenes;
    }

    @NotNull
    public Set<String> reportableGermlineGenes() {
        return reportableGenes;
    }

    @NotNull
    Set<String> notifiableGenes(@NotNull LimsGermlineReportingLevel reportingLevel) {
        if (reportingLevel.equals(LimsGermlineReportingLevel.REPORT_WITH_NOTIFICATION)) {
            return notifiableGenes;
        }

        return Collections.emptySet();
    }

    public boolean notifyAboutGene(@NotNull LimsGermlineReportingLevel reportingLevel, @NotNull String germlineGene) {
        boolean reportableContains = reportableGenes.contains(germlineGene);
        if (!reportableContains) {
            LOGGER.warn("Requested notification status for a gene that is not amongst set of reportable germline genes: {}", germlineGene);
        }

        return reportingLevel.equals(LimsGermlineReportingLevel.REPORT_WITH_NOTIFICATION) &&  notifiableGenes.contains(germlineGene);
    }
}
