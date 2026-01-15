package com.hartwig.hmftools.datamodel.finding;

import java.util.List;

import com.hartwig.hmftools.datamodel.chord.ChordStatus;

import org.jetbrains.annotations.NotNull;

import io.soabase.recordbuilder.core.RecordBuilder;

@RecordBuilder
public record HomologousRecombination(
        @NotNull String findingKey,
        double brca1Value,
        double brca2Value,
        double hrdValue,
        @NotNull ChordStatus hrStatus,
        @NotNull String hrdType,
        @NotNull List<GainDeletion> lohCopyNumbers,
        @NotNull List<String> genes
) implements Finding {}
