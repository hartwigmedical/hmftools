package com.hartwig.hmftools.compar.common;

public enum SourceType
{
    OLD,
    NEW;

    public String configStr() { return this.toString().toLowerCase(); }
}
