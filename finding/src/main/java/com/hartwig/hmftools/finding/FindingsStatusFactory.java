package com.hartwig.hmftools.finding;

import java.util.Objects;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;

import com.hartwig.hmftools.datamodel.purple.PurpleQCStatus;
import com.hartwig.hmftools.finding.datamodel.finding.FindingStatus;
import com.hartwig.hmftools.finding.datamodel.finding.FindingStatusBuilder;

import org.jetbrains.annotations.Nullable;

class FindingsStatusFactory
{
    static FindingStatus toFindingsStatus(Set<PurpleQCStatus> purpleQCStatuses)
    {
        SortedSet<FindingStatus.ResultIssue> errors = convert(errors(purpleQCStatuses));
        return FindingStatusBuilder.builder()
                .status(errors.isEmpty() ? FindingStatus.ResultStatus.OK : FindingStatus.ResultStatus.NOT_RELIABLE)
                .errors(errors)
                .warnings(convert(warnings(purpleQCStatuses)))
                .build();
    }

    private static Set<PurpleQCStatus> warnings(Set<PurpleQCStatus> purpleQCStatuses)
    {
        return filter(purpleQCStatuses, "WARN_");
    }

    private static Set<PurpleQCStatus> errors(Set<PurpleQCStatus> purpleQCStatuses)
    {
        return filter(purpleQCStatuses, "FAIL_");
    }

    private static Set<PurpleQCStatus> filter(Set<PurpleQCStatus> purpleQCStatuses, String prefix)
    {
        return purpleQCStatuses.stream().filter(s -> s.name().startsWith(prefix)).collect(Collectors.toSet());
    }

    private static SortedSet<FindingStatus.ResultIssue> convert(Set<PurpleQCStatus> purpleQCStatuses)
    {
        return purpleQCStatuses.stream()
                .map(FindingsStatusFactory::convert)
                .filter(Objects::nonNull)
                .collect(Collectors.toCollection(TreeSet::new));
    }

    @Nullable
    private static FindingStatus.ResultIssue convert(PurpleQCStatus purpleQCStatus)
    {
        return switch(purpleQCStatus)
        {
            case PASS -> null;
            case WARN_DELETED_GENES -> FindingStatus.ResultIssue.DELETED_GENES;
            case WARN_HIGH_COPY_NUMBER_NOISE -> FindingStatus.ResultIssue.HIGH_COPY_NUMBER_NOISE;
            case WARN_GENDER_MISMATCH -> FindingStatus.ResultIssue.GENDER_MISMATCH;
            case WARN_LOW_PURITY -> FindingStatus.ResultIssue.LOW_PURITY;
            case WARN_TINC, FAIL_TINC -> FindingStatus.ResultIssue.TUMOR_IN_NORMAL_CONTAMINATION;
            case FAIL_CONTAMINATION -> FindingStatus.ResultIssue.CONTAMINATION;
            case FAIL_NO_TUMOR -> FindingStatus.ResultIssue.NO_TUMOR;
        };
    }
}
