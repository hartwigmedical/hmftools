package com.hartwig.hmftools.patientreporter.copynumber;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.actionability.cancertype.CancerTypeMappingReading;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.patientreporter.genepanel.GeneModel;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class CopyNumberAnalyzer {

    private CopyNumberAnalyzer() {
    }

    @NotNull
    public static CopyNumberAnalysis analyzeCopyNumbers(@NotNull PurityContext purityContext,
            @NotNull List<PurpleCopyNumber> purpleCopyNumbers, @NotNull List<GeneCopyNumber> exomeGeneCopyNumbers,
            @NotNull ActionabilityAnalyzer actionabilityAnalyzer, @Nullable PatientTumorLocation patientTumorLocation,
            @NotNull GeneModel geneModel) throws IOException {
        String primaryTumorLocation = patientTumorLocation != null ? patientTumorLocation.primaryTumorLocation() : Strings.EMPTY;
        CancerTypeMappingReading cancerTypeMappingReading = CancerTypeMappingReading.readingFile();
        String doidsPrimaryTumorLocation = cancerTypeMappingReading.doidsForPrimaryTumorLocation(primaryTumorLocation);

        List<GeneCopyNumber> panelGeneCopyNumbers = exomeGeneCopyNumbers.stream()
                .filter(geneCopyNumber -> geneModel.cnvGenePanel().contains(geneCopyNumber.gene()))
                .collect(Collectors.toList());

        FittedPurity bestFit = purityContext.bestFit();
        List<GeneCopyNumber> reportableGeneCopyNumbers =
                ReportableCopyNumbers.filterCopyNumbersForReporting(panelGeneCopyNumbers, bestFit.ploidy(), geneModel);

        Map<GeneCopyNumber, List<EvidenceItem>> evidencePerGeneCopyNumber = findEvidenceForCopyNumbers(exomeGeneCopyNumbers,
                doidsPrimaryTumorLocation,
                actionabilityAnalyzer,
                bestFit.ploidy(),
                geneModel);

        // KODU: Add gene copy numbers for which evidence has been found but which were not selected yet.
        for (GeneCopyNumber actionableGeneCopyNumber : evidencePerGeneCopyNumber.keySet()) {
            if (!evidencePerGeneCopyNumber.get(actionableGeneCopyNumber).isEmpty() && !reportableGeneCopyNumbers.contains(
                    actionableGeneCopyNumber)) {
                reportableGeneCopyNumbers.add(actionableGeneCopyNumber);
            }
        }

        return ImmutableCopyNumberAnalysis.builder()
                .gender(purityContext.gender())
                .status(purityContext.status())
                .fittedPurity(bestFit)
                .fittedScorePurity(purityContext.score())
                .copyNumbers(purpleCopyNumbers)
                .exomeGeneCopyNumbers(exomeGeneCopyNumbers)
                .reportableGeneCopyNumbers(reportableGeneCopyNumbers)
                .evidenceItems(toList(evidencePerGeneCopyNumber))
                .build();
    }

    @NotNull
    private static List<EvidenceItem> toList(@NotNull Map<GeneCopyNumber, List<EvidenceItem>> evidencePerCopyNumber) {
        List<EvidenceItem> evidenceItemList = Lists.newArrayList();
        for (List<EvidenceItem> itemsPerVariant : evidencePerCopyNumber.values()) {
            evidenceItemList.addAll(itemsPerVariant);
        }
        return evidenceItemList;
    }

    @NotNull
    private static Map<GeneCopyNumber, List<EvidenceItem>> findEvidenceForCopyNumbers(@NotNull List<GeneCopyNumber> geneCopyNumbers,
            @Nullable String doidsPrimaryTumorLocation, @NotNull ActionabilityAnalyzer actionabilityAnalyzer, double averageTumorPloidy,
            @NotNull GeneModel geneModel) {
        Set<String> actionableGenesCNVS = actionabilityAnalyzer.cnvAnalyzer().actionableGenes();
        Map<GeneCopyNumber, List<EvidenceItem>> evidenceItemsCopyNumber = Maps.newHashMap();

        List<GeneCopyNumber> geneCopyNumbersOnActionableGenes = geneCopyNumbers.stream()
                .filter(geneCopyNumber -> actionableGenesCNVS.contains(geneCopyNumber.gene()))
                .collect(Collectors.toList());

        List<GeneCopyNumber> reportableGeneCopyNumbers =
                ReportableCopyNumbers.filterCopyNumbersForReporting(geneCopyNumbersOnActionableGenes, averageTumorPloidy, geneModel);

        for (GeneCopyNumber geneCopyNumber : reportableGeneCopyNumbers) {
            evidenceItemsCopyNumber.put(geneCopyNumber,
                    actionabilityAnalyzer.cnvAnalyzer()
                            .evidenceForCopyNumberEvent(geneCopyNumber,
                                    doidsPrimaryTumorLocation,
                                    actionabilityAnalyzer.cancerTypeAnalyzer()));
        }

        return evidenceItemsCopyNumber;
    }
}
