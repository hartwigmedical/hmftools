package com.hartwig.hmftools.sage.common;

import com.hartwig.hmftools.common.variant.VariantReadSupport;

public enum ReadContextMatch
{
    NONE(false, false, false),
    CORE(false, true, false), // matches core only
    PARTIAL_CORE(false, true, true), // matches part of core and all of one flank
    FULL(false, true, true), // matches core and both flanks entirely
    REALIGNED(false, true, true),
    REF(true, false, false), // matches the ref at least for the core
    SIMPLE_ALT(false, false, false); // matches a SNV/MNV without the core matching

    public final boolean SupportsAlt;
    public final boolean SupportsRef;
    public final boolean FullAltSupport;

    ReadContextMatch(boolean supportsRef, boolean supportsAlt, boolean fullAltSupport)
    {
        SupportsRef = supportsRef;
        SupportsAlt = supportsAlt;
        FullAltSupport = fullAltSupport;
    }

    public VariantReadSupport toReadSupport()
    {
        switch(this)
        {
            case FULL: return VariantReadSupport.FULL;
            case PARTIAL_CORE: return VariantReadSupport.PARTIAL_CORE;
            case REALIGNED: return VariantReadSupport.REALIGNED;
            case CORE: return VariantReadSupport.CORE;
            case REF: return VariantReadSupport.REF;
            default: return null;
        }
    }

}
