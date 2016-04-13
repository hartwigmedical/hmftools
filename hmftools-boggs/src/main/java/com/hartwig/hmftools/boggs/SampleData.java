package com.hartwig.hmftools.boggs;

import com.hartwig.hmftools.boggs.flagstatreader.FlagstatData;
import org.jetbrains.annotations.NotNull;

public class SampleData {

    @NotNull
    private final String sampleID = "IRRELEVANT"; // KODU: Eventually should contain sample barcode (FR12345678)
    @NotNull
    private final String externalID;
    @NotNull
    private final Iterable<FlagstatData> mappingFlagstats;
    @NotNull
    private final FlagstatData markdupFlagstatData;
    @NotNull
    private final FlagstatData realignFlagstatData;

    public SampleData(@NotNull String externalID, @NotNull Iterable<FlagstatData> mappingFlagstats,
                      @NotNull FlagstatData markdupFlagstatData, @NotNull FlagstatData realignFlagstatData) {
        this.externalID = externalID;
        this.mappingFlagstats = mappingFlagstats;
        this.markdupFlagstatData = markdupFlagstatData;
        this.realignFlagstatData = realignFlagstatData;
    }

    @NotNull
    public String sampleID() {
        return sampleID;
    }

    @NotNull
    public String externalID() {
        return externalID;
    }

    @NotNull
    public Iterable<FlagstatData> mappingFlagstats() {
        return mappingFlagstats;
    }

    @NotNull
    public FlagstatData markdupFlagstat() {
        return markdupFlagstatData;
    }

    @NotNull
    public FlagstatData realignFlagstat() {
        return realignFlagstatData;
    }

    @Override
    public String toString() {
        return "SampleData{" +
                "sampleID='" + sampleID + '\'' +
                ", externalID='" + externalID + '\'' +
                ", mappingFlagstats=" + mappingFlagstats +
                ", markdupFlagstatData=" + markdupFlagstatData +
                ", realignFlagstatData=" + realignFlagstatData +
                '}';
    }
}
