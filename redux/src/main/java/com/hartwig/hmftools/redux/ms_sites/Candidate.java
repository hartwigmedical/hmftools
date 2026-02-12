package com.hartwig.hmftools.redux.ms_sites;

import htsjdk.samtools.util.StringUtil;

public class Candidate
{
    int startIndex;

    // where we are up to, it is also same as end
    int currentEndIndex;
    byte[] pattern;

    String patternString;

    boolean complete = false;

    Candidate(byte[] pattern, int startIndex, int endIndex)
    {
        this.pattern = pattern;
        this.startIndex = startIndex;
        this.currentEndIndex = endIndex;
        patternString = StringUtil.bytesToString(pattern);
    }

    byte nextBase()
    {
        return pattern[(currentEndIndex - startIndex) % pattern.length];
    }

    int length()
    {
        return currentEndIndex - startIndex;
    }

    int numFullUnits()
    {
        return length() / pattern.length;
    }
}
