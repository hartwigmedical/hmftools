package com.hartwig.hmftools.finding;

import java.util.Objects;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;

import com.hartwig.hmftools.finding.datamodel.Qc;
import com.hartwig.hmftools.finding.datamodel.finding.FindingStatus;
import com.hartwig.hmftools.finding.datamodel.finding.FindingStatusBuilder;

import org.jetbrains.annotations.Nullable;

class FindingsStatusFactory
{
    static FindingStatus toFindingsStatus(Set<Qc.QCStatus> qcStatuses)
    {
        SortedSet<FindingStatus.Issue> errors = convert(errors(qcStatuses));
        return FindingStatusBuilder.builder()
                .status(errors.isEmpty() ? FindingStatus.Status.OK : FindingStatus.Status.NOT_RELIABLE)
                .errors(errors)
                .warnings(convert(warnings(qcStatuses)))
                .build();
    }

    private static Set<Qc.QCStatus> warnings(Set<Qc.QCStatus> qcStatuses)
    {
        return filter(qcStatuses, "WARN_");
    }

    private static Set<Qc.QCStatus> errors(Set<Qc.QCStatus> qcStatuses)
    {
        return filter(qcStatuses, "FAIL_");
    }

    private static Set<Qc.QCStatus> filter(Set<Qc.QCStatus> qcStatuses, String prefix)
    {
        return qcStatuses.stream().filter(s -> s.name().startsWith(prefix)).collect(Collectors.toSet());
    }

    private static SortedSet<FindingStatus.Issue> convert(Set<Qc.QCStatus> qcStatuses)
    {
        return qcStatuses.stream()
                .map(FindingsStatusFactory::convert)
                .filter(Objects::nonNull)
                .collect(Collectors.toCollection(TreeSet::new));
    }

    @Nullable
    private static FindingStatus.Issue convert(Qc.QCStatus qcStatus)
    {
        return switch(qcStatus)
        {
            case PASS -> null;
            case WARN_DELETED_GENES -> FindingStatus.Issue.DELETED_GENES;
            case WARN_HIGH_COPY_NUMBER_NOISE -> FindingStatus.Issue.HIGH_COPY_NUMBER_NOISE;
            case WARN_GENDER_MISMATCH -> FindingStatus.Issue.GENDER_MISMATCH;
            case WARN_LOW_PURITY -> FindingStatus.Issue.LOW_PURITY;
            case WARN_TUMOR_IN_NORMAL_CONTAMINATION, FAIL_TUMOR_IN_NORMAL_CONTAMINATION ->
                    FindingStatus.Issue.TUMOR_IN_NORMAL_CONTAMINATION;
            case WARN_TUMOR_LOW_COVERAGE -> FindingStatus.Issue.TUMOR_SAMPLE_QUALITY_CONTROL;
            case WARN_TUMOR_LOW_MAPPED_PROPORTION -> FindingStatus.Issue.TUMOR_SAMPLE_QUALITY_CONTROL;
            case WARN_TUMOR_LOW_BASE_QUAL -> FindingStatus.Issue.TUMOR_SAMPLE_QUALITY_CONTROL;
            case WARN_TUMOR_LOW_MAP_QUAL -> FindingStatus.Issue.TUMOR_SAMPLE_QUALITY_CONTROL;
            case WARN_NORMAL_LOW_COVERAGE -> FindingStatus.Issue.REF_SAMPLE_QUALITY_CONTROL;
            case WARN_NORMAL_LOW_MAPPED_PROPORTION -> FindingStatus.Issue.REF_SAMPLE_QUALITY_CONTROL;
            case WARN_NORMAL_LOW_BASE_QUAL -> FindingStatus.Issue.REF_SAMPLE_QUALITY_CONTROL;
            case WARN_NORMAL_LOW_MAP_QUAL -> FindingStatus.Issue.REF_SAMPLE_QUALITY_CONTROL;
            case FAIL_CONTAMINATION -> FindingStatus.Issue.CONTAMINATION;
            case FAIL_NO_TUMOR -> FindingStatus.Issue.NO_TUMOR;
        };
    }
}
