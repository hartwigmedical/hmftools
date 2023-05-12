package com.hartwig.hmftools.orange.conversion;

import java.util.List;

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
        AnnotatedVirus virus = minimalAnnotatedVirusBuilder().build();
        return ImmutableVirusInterpreterData.builder().addReportableViruses(virus).addAllAllViruses(List.of(virus)).build();
    }

    static VirusInterpreterData exhaustiveVirusInterpreterData() {
        AnnotatedVirus virus = minimalAnnotatedVirusBuilder().interpretation("MCV").expectedClonalCoverage(1.0).build();
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
