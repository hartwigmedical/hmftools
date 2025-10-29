package com.hartwig.hmftools.compar.linx;

import static com.hartwig.hmftools.compar.common.CommonUtils.FLD_REPORTED;
import static com.hartwig.hmftools.compar.linx.FusionData.FLD_CHAIN_LINKS;
import static com.hartwig.hmftools.compar.linx.FusionData.FLD_CHAIN_TERM;
import static com.hartwig.hmftools.compar.linx.FusionData.FLD_DOMAINS_KEPT;
import static com.hartwig.hmftools.compar.linx.FusionData.FLD_DOMAINS_LOST;
import static com.hartwig.hmftools.compar.linx.FusionData.FLD_EXON_DOWN;
import static com.hartwig.hmftools.compar.linx.FusionData.FLD_EXON_UP;
import static com.hartwig.hmftools.compar.linx.FusionData.FLD_JUNCTION_COPY_NUMBER;
import static com.hartwig.hmftools.compar.linx.FusionData.FLD_LIKELIHOOD;
import static com.hartwig.hmftools.compar.linx.FusionData.FLD_PHASED;
import static com.hartwig.hmftools.compar.linx.FusionData.FLD_REPORTED_TYPE;
import static com.hartwig.hmftools.compar.linx.FusionData.FLD_TRANSCRIPT_DOWN;
import static com.hartwig.hmftools.compar.linx.FusionData.FLD_TRANSCRIPT_UP;

import java.util.HashMap;
import java.util.Map;

import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItemTest;

import org.junit.Before;

public class FusionDataTest extends ComparableItemTest<FusionData, FusionComparer, TestFusionDataBuilder>
{
    @Before
    public void setUp()
    {
        comparer = new FusionComparer(new ComparConfig());
        builder = TestFusionDataBuilder.BUILDER;

        FusionData alternateValueSource = builder.createWithAlternateDefaults();

        fieldToAlternateValueInitializer = new HashMap<>();
        fieldToAlternateValueInitializer.put(FLD_REPORTED_TYPE, b -> b.reportedType = alternateValueSource.Fusion.reportedType());
        fieldToAlternateValueInitializer.put(FLD_PHASED, b -> b.phased = alternateValueSource.Fusion.phased());
        fieldToAlternateValueInitializer.put(FLD_LIKELIHOOD, b -> b.likelihood = alternateValueSource.Fusion.likelihood());
        fieldToAlternateValueInitializer.put(FLD_TRANSCRIPT_UP, b -> b.transcriptUp = alternateValueSource.Fusion.geneTranscriptStart());
        fieldToAlternateValueInitializer.put(FLD_EXON_UP, b -> b.exonUp = alternateValueSource.Fusion.fusedExonUp());
        fieldToAlternateValueInitializer.put(FLD_TRANSCRIPT_DOWN, b -> b.transcriptDown = alternateValueSource.Fusion.geneTranscriptEnd());
        fieldToAlternateValueInitializer.put(FLD_EXON_DOWN, b -> b.exonDown = alternateValueSource.Fusion.fusedExonDown());
        fieldToAlternateValueInitializer.put(FLD_CHAIN_LINKS, b -> b.chainLinks = alternateValueSource.Fusion.chainLinks());
        fieldToAlternateValueInitializer.put(FLD_CHAIN_TERM, b -> b.chainTerminated = alternateValueSource.Fusion.chainTerminated());
        fieldToAlternateValueInitializer.put(FLD_DOMAINS_KEPT, b -> b.domainsKept = alternateValueSource.Fusion.domainsKept());
        fieldToAlternateValueInitializer.put(FLD_DOMAINS_LOST, b -> b.domainsLost = alternateValueSource.Fusion.domainsLost());
        fieldToAlternateValueInitializer.put(FLD_JUNCTION_COPY_NUMBER, b -> b.junctionCopyNumber =
                alternateValueSource.Fusion.junctionCopyNumber());

        nameToAlternateIndexInitializer = Map.of("FusionName", b -> b.fusionName = alternateValueSource.GeneMappedName);
        reportabilityFieldToFalseReportabilityInitializer = Map.of(FLD_REPORTED, b -> b.reported = false);
    }
}
