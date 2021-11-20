package com.hartwig.hmftools.gripss;

public class VcfIdGenerator
{
    private int mVcfId;

    private static final char SUFFIX_FIRST = 'o';
    private static final char SUFFIX_SECOND = 'h';

    public VcfIdGenerator()
    {
        mVcfId = 1;
    }

    private String formatId(final char suffix) { return String.format("vcf%d%c", mVcfId, suffix); }

    private char suffix(boolean isFirst) { return isFirst ? SUFFIX_FIRST : SUFFIX_SECOND; }

    public String next(boolean isFirst)
    {
        ++mVcfId;
        return formatId(suffix(isFirst));
    }

    public String current(boolean isFirst)
    {
        return formatId(suffix(isFirst));
    }
}
