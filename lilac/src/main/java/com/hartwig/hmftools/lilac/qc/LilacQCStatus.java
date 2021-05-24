/*
 * Decompiled with CFR 0.151.
 *
 * Could not load the following classes:
 *  kotlin.Metadata
 */
package com.hartwig.hmftools.lilac.qc;

public enum LilacQCStatus
{
    PASS,
    WARN_UNMATCHED_TYPE,
    WARN_UNMATCHED_SOMATIC_VARIANT,
    WARN_WILDCARD_MATCH,
    WARN_UNMATCHED_HAPLOTYPE,
    WARN_UNMATCHED_AMINO_ACID,
    WARN_UNMATCHED_INDEL;
}