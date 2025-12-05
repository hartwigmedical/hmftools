package com.hartwig.hmftools.compar.linx;

import static java.util.Collections.emptyList;

import java.util.function.Consumer;

import com.hartwig.hmftools.common.linx.FusionLikelihoodType;
import com.hartwig.hmftools.common.linx.FusionPhasedType;
import com.hartwig.hmftools.common.linx.ImmutableLinxFusion;
import com.hartwig.hmftools.common.linx.LinxFusion;
import com.hartwig.hmftools.compar.TestComparableItemBuilder;

public class TestFusionDataBuilder
{
    public String fusionName = "TMPRSS2_ERG";
    public boolean reported = true;
    public String reportedType = "KNOWN_PAIR";
    public FusionPhasedType phased = FusionPhasedType.INFRAME;
    public FusionLikelihoodType likelihood = FusionLikelihoodType.HIGH;
    public String transcriptUp = "ENST00000332149";
    public int exonUp = 2;
    public String transcriptDown = "ENST00000288319";
    public int exonDown = 3;
    public int chainLinks = 0;
    public boolean chainTerminated = false;
    public String domainsKept = "Pointed domain;Ets domain";
    public String domainsLost = "";
    public Double junctionCopyNumber = 1.0;

    private static final Consumer<TestFusionDataBuilder> ALTERNATE_INITIALIZER = b ->
    {
        b.fusionName = "CSMD3_ZNF704";
        b.reported = false;
        b.reportedType = "NONE";
        b.phased =FusionPhasedType.SKIPPED_EXONS;
        b.likelihood = FusionLikelihoodType.NA;
        b.transcriptUp = "ENST00000297405";
        b.exonUp = 12;
        b.transcriptDown = "ENST00000327835";
        b.exonDown = 7;
        b.chainLinks = 2;
        b.chainTerminated = true;
        b.domainsKept = "Villin headpiece";
        b.domainsLost = "Zinc finger";
        b.junctionCopyNumber = 2.0;
    };

    public static final TestComparableItemBuilder<TestFusionDataBuilder, FusionData> BUILDER =
            new TestComparableItemBuilder<>(TestFusionDataBuilder::new, TestFusionDataBuilder::build, ALTERNATE_INITIALIZER);

    private FusionData build()
    {
        final LinxFusion fusion = ImmutableLinxFusion.builder()
                .name(fusionName)
                .reported(reported)
                .reportedType(reportedType)
                .phased(phased)
                .likelihood(likelihood)
                .geneTranscriptStart(transcriptUp)
                .fusedExonUp(exonUp)
                .geneTranscriptEnd(transcriptDown)
                .fusedExonDown(exonDown)
                .chainLinks(chainLinks)
                .chainTerminated(chainTerminated)
                .domainsKept(domainsKept)
                .domainsLost(domainsLost)
                .junctionCopyNumber(junctionCopyNumber)
                .fivePrimeBreakendId(-1)
                .threePrimeBreakendId(-1)
                .fivePrimeVcfId("")
                .threePrimeVcfId("")
                .fivePrimeCoords("")
                .threePrimeCoords("")
                .reportableReasons(emptyList())
                .chainLength(-1)
                .skippedExonsUp(-1)
                .skippedExonsDown(-1)
                .geneStart("")
                .geneContextStart("")
                .geneEnd("")
                .geneContextEnd("")
                .build();
        return new FusionData(fusion, fusionName);
    }
}
