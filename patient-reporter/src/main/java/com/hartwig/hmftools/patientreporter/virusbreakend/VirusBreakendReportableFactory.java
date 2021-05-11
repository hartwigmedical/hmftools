package com.hartwig.hmftools.patientreporter.virusbreakend;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.virusbreakend.VirusBreakend;
import com.hartwig.hmftools.common.virusbreakend.VirusBreakendQCStatus;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class VirusBreakendReportableFactory {

    private VirusBreakendReportableFactory() {

    }

    @NotNull
    public static ReportableVirusBreakendTotal analyzeVirusBreakend(@NotNull List<VirusBreakend> virusBreakends,
            @NotNull VirusDbModel virusDbModel, @NotNull VirusSummaryModel virusSummaryModel,
            @NotNull VirusBlackListModel virusBlackListModel) {

        List<VirusBreakend> virusBreakendsFiltered = Lists.newArrayList();
        List<ReportableVirusbreakend> virusBreakendsReportable = Lists.newArrayList();
        Set<String> postiveSummary = Sets.newHashSet();
        Set<String> negativeSummary = Sets.newHashSet();
        Set<String> summary = Sets.newHashSet();

        for (VirusBreakend virusBreakend : virusBreakends) {
            if (virusBreakend.QCStatus() != VirusBreakendQCStatus.LOW_VIRAL_COVERAGE) {
                if (virusBreakend.integrations() >= 1) {
                    if (virusBlackListModel.checkTaxusForId(virusBreakend.taxidGenus()) != null && virusBlackListModel.checkTaxusForId(
                            virusBreakend.taxidGenus()).equals("taxid_genus")) {
                        if (!virusBlackListModel.checkVirusForBlacklisting(virusBreakend.taxidGenus())) {
                            virusBreakendsFiltered.add(virusBreakend);
                        }
                    } else if (virusBlackListModel.checkTaxusForId(virusBreakend.taxidGenus()) != null
                            && virusBlackListModel.checkTaxusForId(virusBreakend.taxidGenus()).equals("taxid_species")) {
                        if (!virusBlackListModel.checkVirusForBlacklisting(virusBreakend.taxidSpecies())) {
                            virusBreakendsFiltered.add(virusBreakend);
                        }
                    } else {
                        virusBreakendsFiltered.add(virusBreakend);
                    }
                }
            }
        }

        for (VirusBreakend virusBreakend : virusBreakendsFiltered) {

            String virusName = virusDbModel.findVirus(virusBreakend.referenceTaxid());
            virusBreakendsReportable.add(ImmutableReportableVirusbreakend.builder()
                    .virusName(virusName)
                    .integrations(virusBreakend.integrations())
                    .build());

            if (virusSummaryModel.mapIdtoVirusName(virusBreakend.taxidSpecies())) {
                if (!virusSummaryModel.findVirusSummary(virusBreakend.taxidSpecies()).equals(Strings.EMPTY)) {
                    postiveSummary.add(virusSummaryModel.findVirusSummary(virusBreakend.taxidSpecies()) + " positive");
                }
            }
        }

        for (VirusBreakend virusBreakend : virusBreakends) {
            for (String virus : virusSummaryModel.virussen()) {
                if (!postiveSummary.contains(virus + " positive")) {
                    if (!virusSummaryModel.findVirusSummary(virusBreakend.taxidSpecies()).equals(Strings.EMPTY)) {
                        negativeSummary.add(virusSummaryModel.findVirusSummary(virusBreakend.taxidSpecies()) + " negative");
                    }
                }
            }
        }

        summary.addAll(postiveSummary);
        summary.addAll(negativeSummary);
        String virusNameSummary = summary.toString().replace("[", "").replace("]", "");

        return ImmutableReportableVirusBreakendTotal.builder()
                .reportableVirussen(virusBreakendsReportable)
                .virusNameSummary(virusNameSummary)
                .build();
    }
}
