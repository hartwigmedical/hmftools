package com.hartwig.hmftools.breakpointinspector.datamodel;

public class BreakpointStats {

    private int prOnlyNormal = 0;
    private int prSrNormal = 0;
    private int prOnlySupport = 0;
    private int srOnlySupport = 0;
    private int prSrSupport = 0;

    public void incrementPrOnlyNormal() {
        this.prOnlyNormal++;
    }

    public void incrementPrSrNormal() {
        this.prSrNormal++;
    }

    public void incrementPrOnlySupport() {
        this.prOnlySupport++;
    }

    public void incrementSrOnlySupport() {
        this.srOnlySupport++;
    }

    public void incrementPrSrSupport() {
        this.prSrSupport++;
    }

    public int prOnlyNormal() {
        return prOnlyNormal;
    }

    public int prSrNormal() {
        return prSrNormal;
    }

    public int prOnlySupport() {
        return prOnlySupport;
    }

    public int srOnlySupport() {
        return srOnlySupport;
    }

    public int prSrSupport() {
        return prSrSupport;
    }
}
