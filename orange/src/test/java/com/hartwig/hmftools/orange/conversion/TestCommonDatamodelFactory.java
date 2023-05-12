package com.hartwig.hmftools.orange.conversion;

import com.hartwig.hmftools.common.virus.AnnotatedVirus;
import com.hartwig.hmftools.common.virus.ImmutableAnnotatedVirus;
import com.hartwig.hmftools.common.virus.ImmutableVirusInterpreterData;
import com.hartwig.hmftools.common.virus.VirusBreakendQCStatus;
import com.hartwig.hmftools.common.virus.VirusInterpreterData;
import com.hartwig.hmftools.common.virus.VirusLikelihoodType;

class TestCommonDatamodelFactory {

    static VirusInterpreterData emptyVirusInterpreterData() {
        return ImmutableVirusInterpreterData.builder().build();
    }

    static VirusInterpreterData minimalVirusInterpreterData() {
        return virusInterpreterDataWith(minimalAnnotatedVirusBuilder().build());
    }

    static VirusInterpreterData exhaustiveVirusInterpreterData() {
        return virusInterpreterDataWith(minimalAnnotatedVirusBuilder().interpretation("MCV").expectedClonalCoverage(1.0).build());
    }

    private static ImmutableVirusInterpreterData virusInterpreterDataWith(AnnotatedVirus virus) {
        return ImmutableVirusInterpreterData.builder().addReportableViruses(virus).addAllViruses(virus).build();
    }

    private static ImmutableAnnotatedVirus.Builder minimalAnnotatedVirusBuilder() {
        return ImmutableAnnotatedVirus.builder()
                .taxid(1)
                .name("virus_name")
                .qcStatus(VirusBreakendQCStatus.NO_ABNORMALITIES)
                .integrations(1)
                .virusDriverLikelihoodType(VirusLikelihoodType.HIGH)
                .percentageCovered(1.0)
                .meanCoverage(1.0)
                .reported(true);
    }
}
