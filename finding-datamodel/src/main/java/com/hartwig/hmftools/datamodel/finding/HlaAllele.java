package com.hartwig.hmftools.datamodel.finding;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import io.soabase.recordbuilder.core.RecordBuilder;

@RecordBuilder
public record HlaAllele(
        @NotNull String findingKey,
        @NotNull String gene,
        @NotNull String allele,
        int germlineCopyNumber,
        @Nullable Double tumorCopyNumber,
        @Nullable Integer refFragments,
        int tumorFragments,
        @Nullable Integer rnaFragments,
        double somaticMissense,
        double somaticNonsenseOrFrameshift,
        double somaticSplice,
        double somaticSynonymous,
        double somaticInframeIndel
) implements Finding {}
