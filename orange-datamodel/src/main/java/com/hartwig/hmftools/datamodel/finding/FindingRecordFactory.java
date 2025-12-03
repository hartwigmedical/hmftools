package com.hartwig.hmftools.datamodel.finding;

import java.util.Optional;

import com.hartwig.hmftools.datamodel.chord.ChordRecord;
import com.hartwig.hmftools.datamodel.cuppa.CuppaData;
import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterData;

// to reduce duplication, the findings are collected from
// various part of the orange record
public class FindingRecordFactory
{
    public static FindingRecord fromOrangeRecord(OrangeRecord orangeRecord)
    {
        return ImmutableFindingRecord.builder()
                .driverSomaticSmallVariants(orangeRecord.purple().driverSomaticSmallVariants())
                .driverGermlineSmallVariants(orangeRecord.purple().driverGermlineSmallVariants())
                .driverSomaticFusions(orangeRecord.linx().reportableSomaticFusions())
                .driverSomaticDisruptions(orangeRecord.linx().driverSomaticDisruptions())
                .driverGermlineDisruptions(orangeRecord.linx().driverGermlineDisruptions())
                .driverViruses(Optional.ofNullable(orangeRecord.virusInterpreter()).map(VirusInterpreterData::driverViruses).orElse(null))
                .microsatelliteStability(orangeRecord.purple().characteristics().microsatelliteStability())
                .tumorMutationStatus(orangeRecord.purple().characteristics().tumorMutationStatus())
                .predictedTumorOrigin(Optional.ofNullable(orangeRecord.cuppa()).map(CuppaData::bestPrediction).orElse(null))
                .homologousRecombination(Optional.ofNullable(orangeRecord.chord()).map(ChordRecord::homologousRecombination).orElse(null))
                .build();
    }
}
