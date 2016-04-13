package com.hartwig.hmftools.boggs;

import com.hartwig.hmftools.boggs.flagstatreader.FlagStatData;
import org.jetbrains.annotations.NotNull;

public class SampleData {

    @NotNull
    private final String sampleID = "IRRELEVANT"; // KODU: Eventually should contain sample barcode (FR12345678)
    @NotNull
    private final String externalID;
    @NotNull
    private final Iterable<FlagStatData> mappingFlagstats;
    @NotNull
    private final FlagStatData markdupFlagstatData;
    @NotNull
    private final FlagStatData realignFlagstatData;

    public SampleData(@NotNull String externalID, @NotNull Iterable<FlagStatData> mappingFlagstats,
                      @NotNull FlagStatData markdupFlagstatData, @NotNull FlagStatData realignFlagstatData) {
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
    public Iterable<FlagStatData> mappingFlagstats() {
        return mappingFlagstats;
    }

    @NotNull
    public FlagStatData markdupFlagstat() {
        return markdupFlagstatData;
    }

    @NotNull
    public FlagStatData realignFlagstat() {
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
