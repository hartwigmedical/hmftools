package com.hartwig.hmftools.finding.util;

import java.util.Comparator;
import java.util.List;
import java.util.Objects;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.hartwig.hmftools.finding.datamodel.driver.Driver;
import com.hartwig.hmftools.finding.datamodel.driver.DriverFindingList;
import com.hartwig.hmftools.finding.datamodel.driver.DriverFindingListBuilder;
import com.hartwig.hmftools.finding.datamodel.finding.Finding;
import com.hartwig.hmftools.finding.datamodel.finding.FindingList;
import com.hartwig.hmftools.finding.datamodel.finding.FindingListBuilder;
import com.hartwig.hmftools.finding.datamodel.finding.FindingStatus;

import org.jspecify.annotations.Nullable;

import jakarta.validation.constraints.NotNull;

public class FindingsConverter
{
    @NotNull
    public static <I extends Finding, O extends Finding> FindingList<O> convert(@NotNull FindingList<I> findingList,
            Function<FindingStatus, FindingStatus> findingsStatusConverter,
            @Nullable Function<I, O> findingConverter,
            @Nullable Comparator<O> comparator)
    {
        return FindingListBuilder.<O>builder()
                .status(findingsStatusConverter.apply(findingList.status()))
                .findings(convert(findingList.findings(), findingConverter, comparator))
                .build();
    }

    @NotNull
    public static <I extends Driver, O extends Driver> DriverFindingList<O> convert(@NotNull DriverFindingList<I> driverFindingList,
            Function<FindingStatus, FindingStatus> findingsStatusConverter,
            @Nullable Function<I, O> findingConverter,
            @Nullable Comparator<O> comparator)
    {
        return DriverFindingListBuilder.<O>builder()
                .status(findingsStatusConverter.apply(driverFindingList.status()))
                .findings(convert(driverFindingList.findings(), findingConverter, comparator))
                .build();
    }

    static <I, O> List<O> convert(List<I> list, Function<I, O> converter, @Nullable Comparator<O> comparator)
    {
        Stream<O> stream = list.stream().map(converter).filter(Objects::nonNull);
        if(comparator != null)
        {
            stream = stream.sorted(comparator);
        }
        return stream.collect(Collectors.toList());
    }
}
