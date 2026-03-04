package com.hartwig.hmftools.finding.datamodel;

import java.util.List;

import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record FindingList<T extends Finding>(
        @NotNull FindingsStatus status,
        @NotNull List<T> findings) implements IFindingList<T>
{

}
