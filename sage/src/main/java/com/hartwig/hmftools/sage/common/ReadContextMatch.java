package com.hartwig.hmftools.sage.common;

public enum ReadContextMatch
{
    NONE,
    CORE, // matches core only
    PARTIAL, // matches core plus part of a flank
    CORE_PARTIAL, // matches part of core but doesn't overlap it all
    FULL; // matches core and both flanks entirely
}
