package com.hartwig.hmftools.breakpointinspector.datamodel;

import java.util.Arrays;
import java.util.List;

import org.jetbrains.annotations.NotNull;

public class BreakpointStats {

    public int prOnlyNormal = 0;
    public int prSrNormal = 0;
    public int prOnlySupport = 0;
    public int srOnlySupport = 0;
    public int prSrSupport = 0;

    @NotNull
    static List<String> header() {
        return Arrays.asList("PR_ONLY_NORMAL", "PR_SR_NORMAL", "PR_ONLY_SUPPORT", "PR_SR_SUPPORT", "SR_ONLY_SUPPORT");
    }

    @NotNull
    List<Integer> data() {
        return Arrays.asList(prOnlyNormal, prSrNormal, prOnlySupport, prSrSupport, srOnlySupport);
    }
}
