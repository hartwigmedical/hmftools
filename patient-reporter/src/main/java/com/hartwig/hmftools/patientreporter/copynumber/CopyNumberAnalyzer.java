package com.hartwig.hmftools.patientreporter.copynumber;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.drivercatalog.CNADrivers;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.purity.FittedPurityStatus;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.patientreporter.actionability.ReportableEvidenceItemFactory;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class CopyNumberAnalyzer {

    private static final Logger LOGGER = LogManager.getLogger(CopyNumberAnalyzer.class);

    private CopyNumberAnalyzer() {
    }

    @NotNull
    public static CopyNumberAnalysis run(@NotNull PurityContext purityContext, @NotNull List<GeneCopyNumber> exomeGeneCopyNumbers,
            @NotNull ActionabilityAnalyzer actionabilityAnalyzer, @Nullable PatientTumorLocation patientTumorLocation) {
        FittedPurity bestFit = purityContext.bestFit();

        CNADrivers copyNumberDrivers = new CNADrivers();
        List<DriverCatalog> drivers = Lists.newArrayList();
        drivers.addAll(copyNumberDrivers.amplifications(bestFit.ploidy(), exomeGeneCopyNumbers));
        drivers.addAll(copyNumberDrivers.deletions(exomeGeneCopyNumbers));

        String primaryTumorLocation = patientTumorLocation != null ? patientTumorLocation.primaryTumorLocation() : null;
        Map<GeneCopyNumber, List<EvidenceItem>> evidencePerGeneCopyNumber =
                actionabilityAnalyzer.evidenceForCopyNumbers(exomeGeneCopyNumbers, primaryTumorLocation, bestFit.ploidy());

        List<EvidenceItem> filteredEvidence = ReportableEvidenceItemFactory.reportableFlatList(evidencePerGeneCopyNumber);

        List<ReportableGainLoss> reportableGainsAndLosses = toReportableGainsAndLosses(exomeGeneCopyNumbers, drivers);

        // Check that all copy numbers with evidence are reported (since they are in the driver catalog).
        Set<String> reportableGenes = Sets.newHashSet();
        for (ReportableGainLoss gainLoss : reportableGainsAndLosses) {
            reportableGenes.add(gainLoss.gene());
        }

        for (Map.Entry<GeneCopyNumber, List<EvidenceItem>> entry : evidencePerGeneCopyNumber.entrySet()) {
            GeneCopyNumber geneCopyNumber = entry.getKey();
            if (!Collections.disjoint(entry.getValue(), filteredEvidence) && !reportableGenes.contains(geneCopyNumber.gene())) {
                LOGGER.warn("Copy number with evidence not reported: {}!", geneCopyNumber.gene());
            }
        }

        return ImmutableCopyNumberAnalysis.builder()
                .hasReliablePurityFit(purityContext.status() != FittedPurityStatus.NO_TUMOR)
                .purity(bestFit.purity())
                .ploidy(bestFit.ploidy())
                .gender(purityContext.gender())
                .exomeGeneCopyNumbers(exomeGeneCopyNumbers)
                .reportableGainsAndLosses(reportableGainsAndLosses)
                .evidenceItems(filteredEvidence)
                .build();
    }

    @NotNull
    private static List<ReportableGainLoss> toReportableGainsAndLosses(@NotNull List<GeneCopyNumber> exomeGeneCopyNumbers,
            @NotNull List<DriverCatalog> drivers) {
        Set<String> driverGenes = Sets.newHashSet();
        for (DriverCatalog driver : drivers) {
            driverGenes.add(driver.gene());
        }

        Map<String, String> armEndLocusOverrides = buildArmEndLocusOverrides();
        List<ReportableGainLoss> reportableGainsAndLosses = Lists.newArrayList();
        for (GeneCopyNumber copyNumber : exomeGeneCopyNumbers) {
            if (driverGenes.contains(copyNumber.gene())) {
                String armEndLocusOverride = armEndLocusOverrides.get(copyNumber.gene());
                reportableGainsAndLosses.add(ImmutableReportableGainLoss.builder()
                        .chromosome(copyNumber.chromosome())
                        .region(armEndLocusOverride != null ? armEndLocusOverride : copyNumber.chromosomeBand())
                        .gene(armEndLocusOverride != null ? Strings.EMPTY : copyNumber.gene())
                        .interpretation(CopyNumberInterpretation.fromCopyNumber(copyNumber))
                        .copies(Math.round(Math.max(0, copyNumber.minCopyNumber())))
                        .build());
            }
        }

        return reportableGainsAndLosses;
    }

    @NotNull
    private static Map<String, String> buildArmEndLocusOverrides() {
        Map<String, String> armEndLocusOverrides = Maps.newHashMap();
        armEndLocusOverrides.put("AC093642.5", "q telomere");
        armEndLocusOverrides.put("AP001464.4", "q centromere");
        armEndLocusOverrides.put("DOCK8", "p telomere");
        armEndLocusOverrides.put("DUX4L7", "q telomere");
        armEndLocusOverrides.put("LINC01001", "p telomere");
        armEndLocusOverrides.put("OR4A5", "p centromere");
        armEndLocusOverrides.put("OR4F21", "p telomere");
        armEndLocusOverrides.put("PARD6G", "q telomere");
        armEndLocusOverrides.put("PPP2R3B", "p telomere");
        armEndLocusOverrides.put("RP11-417J8.3", "q centromere");
        armEndLocusOverrides.put("SPATA31A7", "q centromere");
        return armEndLocusOverrides;
    }
}
