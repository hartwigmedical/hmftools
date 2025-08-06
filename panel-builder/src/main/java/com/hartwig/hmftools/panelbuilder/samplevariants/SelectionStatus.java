package com.hartwig.hmftools.panelbuilder.samplevariants;

public enum SelectionStatus
{
    NOT_SET,    // Not considered for selection yet
    SELECTED,   // Probe included in panel
    FILTERED,   // Failed internal filters
    // Filters reported to the user
    EXCLUDED_CATEGORY,
    PROXIMATE,
    GENE_LOCATIONS;
}
