package com.hartwig.hmftools.boggs;

import com.hartwig.hmftools.boggs.flagstatreader.Flagstat;
import org.jetbrains.annotations.NotNull;

public class SampleData {

    @NotNull
    private final String sampleID = "IRRELEVANT"; // KODU: Eventually should contain sample barcode (FR12345678)
    @NotNull
    private final String externalID;
    @NotNull
    private final Iterable<Flagstat> mappingFlagstats;
    @NotNull
    private final Flagstat markdupFlagstat;
    @NotNull
    private final Flagstat realignFlagstat;

    public SampleData(@NotNull String externalID, @NotNull Iterable<Flagstat> mappingFlagstats,
                      @NotNull Flagstat markdupFlagstat, @NotNull Flagstat realignFlagstat) {
        this.externalID = externalID;
        this.mappingFlagstats = mappingFlagstats;
        this.markdupFlagstat = markdupFlagstat;
        this.realignFlagstat = realignFlagstat;
    }

    @Override
    public String toString() {
        return "SampleData{" +
                "sampleID='" + sampleID + '\'' +
                ", externalID='" + externalID + '\'' +
                ", mappingFlagstats=" + mappingFlagstats +
                ", markdupFlagstat=" + markdupFlagstat +
                ", realignFlagstat=" + realignFlagstat +
                '}';
    }
}
