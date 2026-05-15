package com.hartwig.hmftools.finding.util;

import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.function.Function;

import com.hartwig.hmftools.finding.datamodel.FindingRecord;
import com.hartwig.hmftools.finding.datamodel.FindingRecordBuilder;
import com.hartwig.hmftools.finding.datamodel.HlaAllele;
import com.hartwig.hmftools.finding.datamodel.HlaAlleleBuilder;
import com.hartwig.hmftools.finding.datamodel.Qc;
import com.hartwig.hmftools.finding.datamodel.QcBuilder;
import com.hartwig.hmftools.finding.datamodel.SequencingScope;
import com.hartwig.hmftools.finding.datamodel.driver.Driver;
import com.hartwig.hmftools.finding.datamodel.driver.DriverFindingList;
import com.hartwig.hmftools.finding.datamodel.driver.DriverFindingListBuilder;
import com.hartwig.hmftools.finding.datamodel.finding.Finding;
import com.hartwig.hmftools.finding.datamodel.finding.FindingItem;
import com.hartwig.hmftools.finding.datamodel.finding.FindingItemBuilder;
import com.hartwig.hmftools.finding.datamodel.finding.FindingList;
import com.hartwig.hmftools.finding.datamodel.finding.FindingListBuilder;
import com.hartwig.hmftools.finding.datamodel.finding.FindingStatus;
import com.hartwig.hmftools.finding.datamodel.finding.FindingStatusBuilder;

import org.jspecify.annotations.Nullable;

import jakarta.validation.constraints.NotNull;

public class LowPurityTransformer
{
    private final static double MIN_PURITY_10_PCT = 0.1;
    private final static double MIN_PURITY_20_PCT = 0.2;
    private final static double MIN_PURITY_30_PCT = 0.3;

    @NotNull
    public static FindingRecord transform(@NotNull FindingRecord record)
    {
        return updateQCLowPurityStatus(updateLowPurityStatus(setLowPurityThresholds(record)));
    }

    private static FindingRecord setLowPurityThresholds(@NotNull FindingRecord record)
    {
        boolean isTargeted = record.metaProperties().sequencingScope() == SequencingScope.TARGETED;
        return FindingRecordBuilder.builder(record)
                .somaticSmallVariants(setPurityThreshold(record.somaticSmallVariants(), isTargeted ? MIN_PURITY_10_PCT : null))
                .somaticDisruptions(setPurityThreshold(record.somaticDisruptions(), MIN_PURITY_20_PCT))
                .somaticGainDeletions(setPurityThreshold(record.somaticGainDeletions(), MIN_PURITY_20_PCT))
                .viruses(setPurityThreshold(record.viruses(), MIN_PURITY_20_PCT))
                .microsatelliteStability(setPurityThreshold(record.microsatelliteStability(), MIN_PURITY_20_PCT))
                .tumorMutationalLoad(setPurityThreshold(record.tumorMutationalLoad(), isTargeted ? MIN_PURITY_10_PCT : MIN_PURITY_20_PCT))
                .tumorMutationalBurden(setPurityThreshold(record.tumorMutationalBurden(), isTargeted
                        ? MIN_PURITY_10_PCT
                        : MIN_PURITY_20_PCT))
                .homologousRecombination(setPurityThreshold(record.homologousRecombination(),
                        isTargeted ? MIN_PURITY_30_PCT : MIN_PURITY_20_PCT))
                .hlaAlleles(setPurityThreshold(record.hlaAlleles(), MIN_PURITY_20_PCT))
                .build();
    }

    private static FindingRecord updateLowPurityStatus(@NotNull FindingRecord record)
    {
        double purity = record.purityPloidyFit().purity();
        return FindingRecordBuilder.builder(record)
                .somaticDisruptions(transform(record.somaticDisruptions(), purity))
                .somaticGainDeletions(transform(record.somaticGainDeletions(), purity))
                .viruses(transform(record.viruses(), purity))
                .microsatelliteStability(transform(record.microsatelliteStability(), purity))
                .tumorMutationalLoad(transform(record.tumorMutationalLoad(), purity))
                .tumorMutationalBurden(transform(record.tumorMutationalBurden(), purity))
                .homologousRecombination(transform(record.homologousRecombination(), purity))
                // For HLA status remains the same, but tumor fields are cleared.
                .hlaAlleles(transform(record.hlaAlleles(), purity, Function.identity(), LowPurityTransformer::transform))
                .build();
    }

    @NotNull
    private static FindingRecord updateQCLowPurityStatus(@NotNull FindingRecord record)
    {
        boolean hasLowPurity = hasLowPurity(record.somaticSmallVariants().status())
                || hasLowPurity(record.somaticGainDeletions().status())
                || hasLowPurity(record.somaticDisruptions().status())
                || hasLowPurity(record.fusions().status())
                || hasLowPurity(record.viruses().status())
                || hasLowPurity(record.chromosomeArmCopyNumbers().status())
                || hasLowPurity(record.hlaAlleles().status())
                || hasLowPurity(record.pharmacoGenotypes().status())
                || hasLowPurity(record.predictedTumorOrigin().status())
                || hasLowPurity(record.microsatelliteStability().status())
                || hasLowPurity(record.tumorMutationalLoad().status())
                || hasLowPurity(record.tumorMutationalBurden().status())
                || hasLowPurity(record.homologousRecombination().status());

        SortedSet<Qc.QCStatus> updatedStatus = new TreeSet<>(record.qc().status());
        if(hasLowPurity)
        {
            updatedStatus.add(Qc.QCStatus.WARN_LOW_PURITY);
        }
        else
        {
            updatedStatus.remove(Qc.QCStatus.WARN_LOW_PURITY);
        }

        Qc updatedQc = QcBuilder.builder(record.qc()).status(updatedStatus).build();
        return FindingRecordBuilder.builder(record).qc(updatedQc).build();
    }

    private static boolean hasLowPurity(@NotNull FindingStatus status)
    {
        return status.errors().contains(FindingStatus.Issue.LOW_PURITY)
                || status.warnings().contains(FindingStatus.Issue.LOW_PURITY);
    }

    @NotNull
    private static <T extends Finding> FindingList<T> setPurityThreshold(@NotNull FindingList<T> findingList,
            @Nullable Double purityThreshold)
    {
        return FindingListBuilder.builder(findingList).purityThreshold(purityThreshold).build();
    }

    @NotNull
    private static <T extends Driver> DriverFindingList<T> setPurityThreshold(@NotNull DriverFindingList<T> driverFindingList,
            @Nullable Double purityThreshold)
    {
        return DriverFindingListBuilder.builder(driverFindingList).purityThreshold(purityThreshold).build();
    }

    @NotNull
    private static <T> FindingItem<T> setPurityThreshold(@NotNull FindingItem<T> findingItem, @Nullable Double purity)
    {
        return FindingItemBuilder.builder(findingItem).purityThreshold(purity).build();
    }

    @NotNull
    private static <T extends Finding> FindingList<T> transform(@NotNull FindingList<T> findingList, double purity,
            @NotNull Function<FindingStatus, FindingStatus> findingsStatusConverter,
            @NotNull Function<T, T> findingConverter)
    {
        if(shouldConvert(findingList.status(), hasLowPurity(purity, findingList.purityThreshold())))
        {
            return FindingRecordTransformerUtil.transformFindingList(findingList, findingsStatusConverter, findingConverter, null);
        }
        else
        {
            return findingList;
        }
    }

    @NotNull
    private static <T extends Driver> DriverFindingList<T> transform(@NotNull DriverFindingList<T> driverFindingList, double purity)
    {
        if(shouldConvert(driverFindingList.status(), hasLowPurity(purity, driverFindingList.purityThreshold())))
        {
            return FindingRecordTransformerUtil.transformDriverFindingList(driverFindingList,
                    LowPurityTransformer::transform,
                    f -> null,
                    null);
        }
        else
        {
            return driverFindingList;
        }
    }

    @NotNull
    private static <T> FindingItem<T> transform(@NotNull FindingItem<T> findingItem, double purity)
    {
        if(shouldConvert(findingItem.status(), hasLowPurity(purity, findingItem.purityThreshold())))
        {
            return FindingItemBuilder.<T>builder()
                    .status(transform(findingItem.status()))
                    .purityThreshold(findingItem.purityThreshold())
                    .build();
        }
        else
        {
            return findingItem;
        }
    }

    private static FindingStatus transform(FindingStatus findingStatus)
    {
        return FindingStatusBuilder.builder()
                .status(FindingStatus.Status.NOT_RELIABLE)
                .errors(addLowPurity(findingStatus.errors()))
                .warnings(removeLowPurity(findingStatus.warnings()))
                .build();
    }

    private static boolean shouldConvert(FindingStatus findingStatus, boolean purity)
    {
        return findingStatus.status() == FindingStatus.Status.OK && purity;
    }

    private static SortedSet<FindingStatus.Issue> addLowPurity(SortedSet<FindingStatus.Issue> issues)
    {
        return FindingUtil.addIssues(issues, Set.of(FindingStatus.Issue.LOW_PURITY));
    }

    private static SortedSet<FindingStatus.Issue> removeLowPurity(SortedSet<FindingStatus.Issue> issues)
    {
        return FindingUtil.removeIssues(issues, Set.of(FindingStatus.Issue.LOW_PURITY));
    }

    private static HlaAllele transform(HlaAllele hlaAllele)
    {
        return HlaAlleleBuilder.builder(hlaAllele)
                .tumorCopyNumber(null)
                .build();
    }

    private static boolean hasLowPurity(double purity, @Nullable Double purityThreshold)
    {
        return purityThreshold != null && purity < purityThreshold;
    }
}
