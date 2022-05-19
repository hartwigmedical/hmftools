package com.hartwig.hmftools.patientreporter.remarks;

import java.util.Map;

import com.google.common.annotations.VisibleForTesting;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class SpecialRemarkModel {

    @NotNull
    private final Map<String, String> sampleToSpecialRemarkMap;

    SpecialRemarkModel(@NotNull final Map<String, String> sampleToSpecialRemarkMap) {
        this.sampleToSpecialRemarkMap = sampleToSpecialRemarkMap;
    }

    @NotNull
    public String findSpecialRemarkForSample(@NotNull String sample) {
        boolean sampleHasSpecialRemark = samplePresentInSpecialRemarks(sample);

        return sampleHasSpecialRemark ? sampleToSpecialRemarkMap.get(sample) : Strings.EMPTY;
    }

    @VisibleForTesting
    int specialRemarkCount() {
        return sampleToSpecialRemarkMap.keySet().size();
    }

    @VisibleForTesting
    boolean samplePresentInSpecialRemarks(@NotNull String sample) {
        return sampleToSpecialRemarkMap.containsKey(sample);
    }
}