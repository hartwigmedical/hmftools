package com.hartwig.hmftools.panelbuilder.samplevariants;

public enum SelectionStatus
{
    NOT_SET,    // Not considered for selection yet
    SELECTED,
    FILTERED,   // Failed internal filters
    EXCEEDS_COUNT,
    PROXIMATE,
    GENE_LOCATIONS;
}
