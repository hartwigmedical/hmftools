package com.hartwig.hmftools.datamodel.purple;

import org.jetbrains.annotations.NotNull;

import java.util.Arrays;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

public enum PurpleVariantEffect {
    STOP_GAINED,
    STOP_LOST,
    START_LOST,
    FRAMESHIFT,
    SPLICE_ACCEPTOR,
    SPLICE_DONOR,
    INFRAME_INSERTION,
    INFRAME_DELETION,
    MISSENSE,
    PHASED_INFRAME_INSERTION,
    PHASED_INFRAME_DELETION,
    SYNONYMOUS,
    INTRONIC,
    FIVE_PRIME_UTR,
    THREE_PRIME_UTR,
    UPSTREAM_GENE,
    NON_CODING_TRANSCRIPT,
    OTHER
}
