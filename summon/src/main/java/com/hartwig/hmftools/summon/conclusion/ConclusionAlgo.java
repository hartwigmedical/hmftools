package com.hartwig.hmftools.summon.conclusion;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.chord.ChordStatus;
import com.hartwig.hmftools.common.linx.ReportableHomozygousDisruption;
import com.hartwig.hmftools.common.purple.copynumber.ReportableGainLoss;
import com.hartwig.hmftools.common.sv.linx.LinxFusion;
import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;
import com.hartwig.hmftools.common.variant.tml.TumorMutationalStatus;
import com.hartwig.hmftools.common.virus.AnnotatedVirus;
import com.hartwig.hmftools.summon.SummonData;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class ConclusionAlgo {

    @NotNull
    public static ActionabilityConclusion generateConclusion(@NotNull SummonData summonData) throws IOException {
        String sampleId = summonData.sampleId();

        List<ReportableVariant> reportableSomaticVariants = summonData.purple().reportableSomaticVariants();
        List<ReportableVariant> reportableGermlineVariants = summonData.purple().reportableGermlineVariants();
        List<ReportableGainLoss> reportableGainLosses = summonData.purple().reportableGainsLosses();
        List<LinxFusion> reportableFusions = summonData.linx().reportableFusions();
        List<ReportableHomozygousDisruption> homozygousDisruptions = summonData.linx().homozygousDisruptions();
        List<AnnotatedVirus> reportableViruses = summonData.virusInterpreter().reportableViruses();

        ChordStatus chordStatus = summonData.chord().hrStatus();
        MicrosatelliteStatus microsatelliteStatus = summonData.purple().microsatelliteStatus();
        TumorMutationalStatus tumorMutationalStatus = summonData.purple().tumorMutationalLoadStatus();
        Double tumorMutationalBurden = summonData.purple().tumorMutationalBurdenPerMb(); //should be > 10

        return ImmutableActionabilityConclusion.builder().conclusion(Strings.EMPTY).build();
    }
}
