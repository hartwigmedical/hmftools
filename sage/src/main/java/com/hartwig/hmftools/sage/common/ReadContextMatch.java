package com.hartwig.hmftools.sage.common;

public enum ReadContextMatch
{
    NONE,
    CORE, // matches core only
    PARTIAL_CORE, // matches part of core and all of one flank
    FULL, // matches core and both flanks entirely
    REF, // matches the ref at least for the core

    // Deprecated
    CORE_PARTIAL,
    PARTIAL;
}
