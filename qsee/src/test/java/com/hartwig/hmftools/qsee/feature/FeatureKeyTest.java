package com.hartwig.hmftools.qsee.feature;

import static org.junit.Assert.assertArrayEquals;

import java.util.List;

import org.junit.Test;

public class FeatureKeyTest
{
    @Test
    public void categoricalFeatureNamesSortedCorrectly()
    {
        List<FeatureKey> testKeys = List.of(
                new FeatureKey("PURITY", FeatureType.SUMMARY_TABLE, SourceTool.PURPLE),

                new FeatureKey("COVERAGE_ABOVE_100", FeatureType.SUMMARY_TABLE, SourceTool.BAM_METRICS),
                new FeatureKey("COVERAGE_ABOVE_20", FeatureType.SUMMARY_TABLE, SourceTool.BAM_METRICS),

                new FeatureKey("DiscordantFragType:INV_SHORT", FeatureType.DISCORDANT_FRAG_FREQ, SourceTool.ESVEE),
                new FeatureKey("DiscordantFragType:DEL_DUP_MEDIUM", FeatureType.DISCORDANT_FRAG_FREQ, SourceTool.ESVEE),
                new FeatureKey("DiscordantFragType:DEL_DUP_SHORT", FeatureType.DISCORDANT_FRAG_FREQ, SourceTool.ESVEE),

                new FeatureKey("ReadCount:>=100", FeatureType.DUPLICATE_FREQ, SourceTool.REDUX),
                new FeatureKey("ReadCount:10", FeatureType.DUPLICATE_FREQ, SourceTool.REDUX),
                new FeatureKey("ReadCount:1", FeatureType.DUPLICATE_FREQ, SourceTool.REDUX),

                new FeatureKey("ConsensusType:NONE;RepeatUnitType:>=3bp repeat;RefNumUnits:4", FeatureType.MS_INDEL_ERROR_RATES, SourceTool.REDUX),
                new FeatureKey("ConsensusType:NONE;RepeatUnitType:2bp repeat;RefNumUnits:4", FeatureType.MS_INDEL_ERROR_RATES, SourceTool.REDUX),
                new FeatureKey("ConsensusType:NONE;RepeatUnitType:C/G repeat;RefNumUnits:4", FeatureType.MS_INDEL_ERROR_RATES, SourceTool.REDUX),
                new FeatureKey("ConsensusType:NONE;RepeatUnitType:A/T repeat;RefNumUnits:4", FeatureType.MS_INDEL_ERROR_RATES, SourceTool.REDUX),

                new FeatureKey("ReadType:NONE;StandardMutation:C>A;OriginalQualBin:HIGH (30+)", FeatureType.BQR_PER_ORIG_QUAL, SourceTool.REDUX),
                new FeatureKey("ReadType:NONE;StandardMutation:C>A;OriginalQualBin:LOW (0-29)", FeatureType.BQR_PER_ORIG_QUAL, SourceTool.REDUX)
        );

        List<String> actualNameOrder = testKeys.stream().sorted().map(FeatureKey::name).toList();

        List<String> expectedNameOrder = List.of(
                "COVERAGE_ABOVE_20",
                "COVERAGE_ABOVE_100",
                "PURITY",

                "DiscordantFragType:DEL_DUP_SHORT",
                "DiscordantFragType:DEL_DUP_MEDIUM",
                "DiscordantFragType:INV_SHORT",

                "ReadCount:1",
                "ReadCount:10",
                "ReadCount:>=100",

                "ReadType:NONE;StandardMutation:C>A;OriginalQualBin:LOW (0-29)",
                "ReadType:NONE;StandardMutation:C>A;OriginalQualBin:HIGH (30+)",

                "ConsensusType:NONE;RepeatUnitType:A/T repeat;RefNumUnits:4",
                "ConsensusType:NONE;RepeatUnitType:C/G repeat;RefNumUnits:4",
                "ConsensusType:NONE;RepeatUnitType:2bp repeat;RefNumUnits:4",
                "ConsensusType:NONE;RepeatUnitType:>=3bp repeat;RefNumUnits:4"
        );

        assertArrayEquals(expectedNameOrder.toArray(), actualNameOrder.toArray());
    }

}
