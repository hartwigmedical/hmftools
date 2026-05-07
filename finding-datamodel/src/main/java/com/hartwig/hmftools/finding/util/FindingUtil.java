package com.hartwig.hmftools.finding.util;

import static com.hartwig.hmftools.finding.datamodel.finding.FindingStatus.Issue.NO_REPORTABLE_VALUE;

import java.util.List;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

import com.hartwig.hmftools.finding.datamodel.driver.Driver;
import com.hartwig.hmftools.finding.datamodel.driver.DriverFindingList;
import com.hartwig.hmftools.finding.datamodel.driver.DriverFindingListBuilder;
import com.hartwig.hmftools.finding.datamodel.finding.FindingItem;
import com.hartwig.hmftools.finding.datamodel.finding.FindingItemBuilder;
import com.hartwig.hmftools.finding.datamodel.finding.FindingStatus;
import com.hartwig.hmftools.finding.datamodel.finding.FindingStatusBuilder;

public class FindingUtil
{
    private static final Set<FindingStatus.Issue> GERMLINE_ISSUES =
            Set.of(FindingStatus.Issue.TUMOR_IN_NORMAL_CONTAMINATION, FindingStatus.Issue.NORMAL_REQUIRED, FindingStatus.Issue.NORMAL_SAMPLE_QUALITY_CONTROL);

    public static FindingStatus somaticStatus(FindingStatus status)
    {
        SortedSet<FindingStatus.Issue> errors = removeIssues(status.errors(), GERMLINE_ISSUES);
        return FindingStatusBuilder.builder()
                .status(errors.isEmpty() ? FindingStatus.Status.OK : status.status())
                .errors(errors)
                .warnings(removeIssues(status.warnings(), GERMLINE_ISSUES))
                .build();
    }

    public static FindingStatus germlineStatus(FindingStatus status)
    {
        SortedSet<FindingStatus.Issue> errors = retainIssues(status.errors(), GERMLINE_ISSUES);
        return FindingStatusBuilder.builder()
                .status(errors.isEmpty() ? FindingStatus.Status.OK : FindingStatus.Status.NOT_RELIABLE)
                .errors(errors)
                .warnings(retainIssues(status.warnings(), GERMLINE_ISSUES))
                .build();
    }

    private static SortedSet<FindingStatus.Issue> retainIssues(Set<FindingStatus.Issue> issues, Set<FindingStatus.Issue> issuesToRetain)
    {
        SortedSet<FindingStatus.Issue> newIssues = new TreeSet<>(issues);
        newIssues.retainAll(issuesToRetain);
        return newIssues;
    }

    public static SortedSet<FindingStatus.Issue> removeIssues(Set<FindingStatus.Issue> issues, Set<FindingStatus.Issue> issuesToRemove)
    {
        SortedSet<FindingStatus.Issue> newIssues = new TreeSet<>(issues);
        newIssues.removeAll(issuesToRemove);
        return newIssues;
    }

    public static SortedSet<FindingStatus.Issue> addIssues(Set<FindingStatus.Issue> issues, Set<FindingStatus.Issue> issuesToRemove)
    {
        SortedSet<FindingStatus.Issue> newIssues = new TreeSet<>(issues);
        newIssues.addAll(issuesToRemove);
        return newIssues;
    }


    public static <T extends Driver> DriverFindingList<T> normalRequired()
    {
        return notAvailableDriverFindingList(Set.of(FindingStatus.Issue.NORMAL_REQUIRED));
    }

    public static <T extends Driver> DriverFindingList<T> notAvailableDriverFindingList(Set<FindingStatus.Issue> errors)
    {
        return DriverFindingListBuilder.<T>builder()
                .status(notAvailableStatus(errors))
                .findings(List.of())
                .build();
    }

    public static <T> FindingItem<T> notAvailableFindingItem(Set<FindingStatus.Issue> errors)
    {
        return FindingItemBuilder.<T>builder()
                .status(notAvailableStatus(errors))
                .build();
    }

    public static FindingStatus notAvailableStatus(Set<FindingStatus.Issue> errors)
    {
        return FindingStatusBuilder.builder()
                .status(FindingStatus.Status.NOT_AVAILABLE)
                .errors(new TreeSet<>(errors))
                .warnings(new TreeSet<>())
                .build();
    }

    public static FindingStatus noReportableValueStatus(FindingStatus findingStatus)
    {
        return noReportableValueStatus(findingStatus, FindingStatus.Status.NOT_AVAILABLE);
    }

    public static FindingStatus noReportableValueStatus(FindingStatus findingStatus, FindingStatus.Status status)
    {
        return FindingStatusBuilder.builder(findingStatus)
                .status(findingStatus.isOK() ? status : findingStatus.status()  )
                .errors(new TreeSet<>(Set.of(NO_REPORTABLE_VALUE)))
                .build();
    }

}
