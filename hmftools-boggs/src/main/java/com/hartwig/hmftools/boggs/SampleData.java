package com.hartwig.hmftools.boggs;

import com.hartwig.hmftools.boggs.flagstatreader.FlagStatData;
import org.jetbrains.annotations.NotNull;

public class SampleData {

    @NotNull
    private final String sampleID = "IRRELEVANT"; // KODU: Eventually should contain sample barcode (FR12345678)
    @NotNull
    private final String externalID;
    @NotNull
    private final Iterable<FlagStatData> rawMappingFlagstats;
    @NotNull
    private final Iterable<FlagStatData> sortedMappingFlagstats;
    @NotNull
    private final FlagStatData markdupFlagstats;
    @NotNull
    private final FlagStatData realignFlagstats;

    public SampleData(@NotNull String externalID, @NotNull Iterable<FlagStatData> rawMappingFlagstats,
                      @NotNull Iterable<FlagStatData> sortedMappingFlagstats, @NotNull FlagStatData markdupFlagstats,
                      @NotNull FlagStatData realignFlagstats) {
        this.externalID = externalID;
        this.rawMappingFlagstats = rawMappingFlagstats;
        this.sortedMappingFlagstats = sortedMappingFlagstats;
        this.markdupFlagstats = markdupFlagstats;
        this.realignFlagstats = realignFlagstats;
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
    public Iterable<FlagStatData> rawMappingFlagstats() {
        return rawMappingFlagstats;
    }

    @NotNull
    public Iterable<FlagStatData> sortedMappingFlagstats() {
        return sortedMappingFlagstats;
    }

    @NotNull
    public FlagStatData markdupFlagstat() {
        return markdupFlagstats;
    }

    @NotNull
    public FlagStatData realignFlagstat() {
        return realignFlagstats;
    }

    @Override
    public String toString() {
        return "SampleData{" +
                "sampleID='" + sampleID + '\'' +
                ", externalID='" + externalID + '\'' +
                ", rawMappingFlagstats=" + rawMappingFlagstats +
                ", sortedMappingFlagstats=" + sortedMappingFlagstats +
                ", markdupFlagstats=" + markdupFlagstats +
                ", realignFlagstats=" + realignFlagstats +
                '}';
    }
}
