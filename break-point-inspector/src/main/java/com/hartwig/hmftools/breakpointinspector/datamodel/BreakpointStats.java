package com.hartwig.hmftools.breakpointinspector.datamodel;

import java.util.Arrays;
import java.util.List;

import org.jetbrains.annotations.NotNull;

public class BreakpointStats {

    public int PR_Only_Normal = 0;
    public int PR_SR_Normal = 0;
    public int PR_Only_Support = 0;
    public int PR_SR_Support = 0;
    public int SR_Only_Support = 0;

    @NotNull
    static List<String> GetHeader() {
        return Arrays.asList("PR_ONLY_NORMAL", "PR_SR_NORMAL", "PR_ONLY_SUPPORT", "PR_SR_SUPPORT", "SR_ONLY_SUPPORT");
    }

    @NotNull
    List<Integer> GetData() {
        return Arrays.asList(PR_Only_Normal, PR_SR_Normal, PR_Only_Support, PR_SR_Support, SR_Only_Support);
    }
}
