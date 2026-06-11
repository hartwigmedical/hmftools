package com.hartwig.hmftools.finding;

import java.util.Arrays;
import java.util.Collection;
import java.util.Objects;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;

import com.hartwig.hmftools.datamodel.purple.PurpleQCStatus;
import com.hartwig.hmftools.finding.datamodel.Qc;

public class QcStatusFactory
{
    private final static Set<PurpleQCStatus> QC_ERROR_STATUSES = errors(Arrays.asList(PurpleQCStatus.values()));
    private final static Set<PurpleQCStatus> QC_WARNING_STATUSES = warnings(Arrays.asList(PurpleQCStatus.values()));

    static SortedSet<Qc.QcStatus> toErrors(Set<PurpleQCStatus> qcStatuses) {
        return QcStatusFactory.toQcStatuses(qcStatuses, QcStatusFactory.QC_ERROR_STATUSES);
    }

    static SortedSet<Qc.QcStatus> toWarnings(Set<PurpleQCStatus> qcStatuses) {
        return QcStatusFactory.toQcStatuses(qcStatuses, QcStatusFactory.QC_WARNING_STATUSES);
    }

    private static SortedSet<Qc.QcStatus> toQcStatuses(Set<PurpleQCStatus> qcStatuses, Set<PurpleQCStatus> statusesToInclude)
    {
        return qcStatuses.stream().filter(statusesToInclude::contains)
                .map(o -> switch(o)
                {
                    case PASS -> null;
                    case WARN_DELETED_GENES -> Qc.QcStatus.DELETED_GENES;
                    case WARN_HIGH_COPY_NUMBER_NOISE -> Qc.QcStatus.HIGH_COPY_NUMBER_NOISE;
                    case WARN_GENDER_MISMATCH -> Qc.QcStatus.GENDER_MISMATCH;
                    case WARN_LOW_PURITY -> Qc.QcStatus.LOW_PURITY;
                    case WARN_TINC, FAIL_TINC -> Qc.QcStatus.TUMOR_IN_NORMAL_CONTAMINATION;
                    case FAIL_CONTAMINATION -> Qc.QcStatus.CONTAMINATION;
                    case FAIL_NO_TUMOR -> Qc.QcStatus.NO_TUMOR;
                })
                .filter(Objects::nonNull).collect(Collectors.toCollection(TreeSet::new));
    }

    private static Set<PurpleQCStatus> warnings(Collection<PurpleQCStatus> qcStatuses)
    {
        return filter(qcStatuses, "WARN_");
    }

    private static Set<PurpleQCStatus> errors(Collection<PurpleQCStatus> qcStatuses)
    {
        return filter(qcStatuses, "FAIL_");
    }

    private static Set<PurpleQCStatus> filter(Collection<PurpleQCStatus> qcStatuses, String prefix)
    {
        return qcStatuses.stream().filter(s -> s.name().startsWith(prefix)).collect(Collectors.toSet());
    }
}
