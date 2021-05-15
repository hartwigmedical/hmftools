package com.hartwig.hmftools.patientreporter.virusbreakend;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.virusbreakend.VirusBreakend;
import com.hartwig.hmftools.common.virusbreakend.VirusBreakendFactory;
import com.hartwig.hmftools.common.virusbreakend.VirusBreakendQCStatus;
import com.hartwig.hmftools.patientreporter.algo.GenomicAnalyzer;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class VirusBreakendReportableFactory {

    private static final Logger LOGGER = LogManager.getLogger(VirusBreakendReportableFactory.class);

    private VirusBreakendReportableFactory() {
    }

    @NotNull
    public static ReportableVirusBreakendTotal analyzeVirusBreakend(@NotNull String virusBreakendTsv,
            @NotNull List<VirusBreakend> virusBreakends, @NotNull VirusDbModel virusDbModel, @NotNull VirusSummaryModel virusSummaryModel,
            @NotNull VirusBlacklistModel virusBlackListModel) {

        List<VirusBreakend> virusBreakendsFiltered = Lists.newArrayList();
        List<ReportableVirusBreakend> virusBreakendsReportable = Lists.newArrayList();
        Set<String> positiveSummary = Sets.newHashSet();
        Set<String> negativeSummary = Sets.newHashSet();
        Set<String> summary = Sets.newHashSet();

        LOGGER.info("Loading virus breakend data from {}", new File(virusBreakendTsv).getParent());
        for (VirusBreakend virusBreakend : virusBreakends) {
            if (virusBreakend.QCStatus() != VirusBreakendQCStatus.LOW_VIRAL_COVERAGE) {
                if (virusBreakend.integrations() >= 1) {
                    if (virusBlackListModel.checkTaxusForId(virusBreakend.taxidGenus()) != null && virusBlackListModel.checkTaxusForId(
                            virusBreakend.taxidGenus()).equals("taxid_genus")) {
                        if (!virusBlackListModel.checkVirusForBlacklisting(virusBreakend.taxidGenus())) {
                            virusBreakendsFiltered.add(virusBreakend);
                        }
                    } else if (virusBlackListModel.checkTaxusForId(virusBreakend.taxidSpecies()) != null
                            && virusBlackListModel.checkTaxusForId(virusBreakend.taxidSpecies()).equals("taxid_species")) {
                        if (!virusBlackListModel.checkVirusForBlacklisting(virusBreakend.taxidSpecies())) {
                            virusBreakendsFiltered.add(virusBreakend);
                        }
                    } else {
                        virusBreakendsFiltered.add(virusBreakend);
                    }
                }
            }
        }

        LOGGER.info(" Loaded {} reportable virus breakend from {}", virusBreakendsFiltered.size(), virusBreakendTsv);

        for (VirusBreakend virusBreakend : virusBreakendsFiltered) {
            String virusName = virusDbModel.findVirus(virusBreakend.referenceTaxid());
            virusBreakendsReportable.add(ImmutableReportableVirusBreakend.builder()
                    .virusName(virusName)
                    .integrations(virusBreakend.integrations())
                    .build());

            if (virusSummaryModel.mapIdToVirusName(virusBreakend.taxidSpecies())) {
                if (!virusSummaryModel.findVirusSummary(virusBreakend.taxidSpecies()).equals(Strings.EMPTY)) {
                    positiveSummary.add(virusSummaryModel.findVirusSummary(virusBreakend.taxidSpecies()) + " positive");
                }
                LOGGER.warn("Virus breakend has called not a HPV/EBV/MCV virus");
            }
        }

        for (String virus : virusSummaryModel.viruses()) {
            if (!positiveSummary.contains(virus + " positive")) {
                negativeSummary.add(virus + " negative");
            }
        }

        summary.addAll(positiveSummary);
        summary.addAll(negativeSummary);
        String virusNameSummary = summary.toString().replace("[", "").replace("]", "");

        return ImmutableReportableVirusBreakendTotal.builder()
                .reportableViruses(virusBreakendsReportable)
                .virusNameSummary(virusNameSummary)
                .build();
    }
}
