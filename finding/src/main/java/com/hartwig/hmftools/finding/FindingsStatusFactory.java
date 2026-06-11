package com.hartwig.hmftools.finding;

import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;

import com.hartwig.hmftools.finding.datamodel.Qc;
import com.hartwig.hmftools.finding.datamodel.finding.FindingStatus;
import com.hartwig.hmftools.finding.datamodel.finding.FindingStatusBuilder;

class FindingsStatusFactory
{
    static FindingStatus toFindingsStatus(Set<Qc.QcStatus> errors, Set<Qc.QcStatus> warnings)
    {
        SortedSet<FindingStatus.Issue> findingStatusErrors = convert(errors);
        return FindingStatusBuilder.builder()
                .status(findingStatusErrors.isEmpty() ? FindingStatus.Status.OK : FindingStatus.Status.NOT_RELIABLE)
                .errors(findingStatusErrors)
                .warnings(convert(warnings))
                .build();
    }

    private static SortedSet<FindingStatus.Issue> convert(Set<Qc.QcStatus> qcStatuses)
    {
        return qcStatuses.stream()
                .map(FindingsStatusFactory::convert)
                .collect(Collectors.toCollection(TreeSet::new));
    }

    private static FindingStatus.Issue convert(Qc.QcStatus qcStatus)
    {
        return switch(qcStatus)
        {
            case DELETED_GENES -> FindingStatus.Issue.DELETED_GENES;
            case HIGH_COPY_NUMBER_NOISE -> FindingStatus.Issue.HIGH_COPY_NUMBER_NOISE;
            case GENDER_MISMATCH -> FindingStatus.Issue.GENDER_MISMATCH;
            case LOW_PURITY -> FindingStatus.Issue.LOW_PURITY;
            case TUMOR_IN_NORMAL_CONTAMINATION -> FindingStatus.Issue.TUMOR_IN_NORMAL_CONTAMINATION;
            case TUMOR_LOW_COVERAGE -> FindingStatus.Issue.TUMOR_SAMPLE_QUALITY_CONTROL;
            case TUMOR_LOW_MAPPED_PROPORTION -> FindingStatus.Issue.TUMOR_SAMPLE_QUALITY_CONTROL;
            case TUMOR_LOW_BASE_QUAL -> FindingStatus.Issue.TUMOR_SAMPLE_QUALITY_CONTROL;
            case TUMOR_LOW_MAP_QUAL -> FindingStatus.Issue.TUMOR_SAMPLE_QUALITY_CONTROL;
            case NORMAL_LOW_COVERAGE -> FindingStatus.Issue.NORMAL_SAMPLE_QUALITY_CONTROL;
            case NORMAL_LOW_MAPPED_PROPORTION -> FindingStatus.Issue.NORMAL_SAMPLE_QUALITY_CONTROL;
            case NORMAL_LOW_BASE_QUAL -> FindingStatus.Issue.NORMAL_SAMPLE_QUALITY_CONTROL;
            case NORMAL_LOW_MAP_QUAL -> FindingStatus.Issue.NORMAL_SAMPLE_QUALITY_CONTROL;
            case CONTAMINATION -> FindingStatus.Issue.CONTAMINATION;
            case NO_TUMOR -> FindingStatus.Issue.NO_TUMOR;
        };
    }
}
