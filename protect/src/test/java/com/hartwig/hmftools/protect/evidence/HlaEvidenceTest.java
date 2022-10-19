package com.hartwig.hmftools.protect.evidence;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.hla.ImmutableLilacSummaryData;
import com.hartwig.hmftools.common.hla.LilacAllele;
import com.hartwig.hmftools.common.hla.LilacSummaryData;
import com.hartwig.hmftools.common.lilac.LilacTestFactory;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.serve.ServeTestFactory;
import com.hartwig.hmftools.common.serve.actionability.immuno.ActionableHLA;
import com.hartwig.hmftools.common.serve.actionability.immuno.ImmutableActionableHLA;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class HlaEvidenceTest {

    @Test
    public void canDetermineEvidenceForHLA() {

        ActionableHLA hla = ImmutableActionableHLA.builder()
                .from(ServeTestFactory.createTestActionableHLA())
                .hlaType("Allele 1")
                .build();

        HlaEvidence hlaEvidence =
                new HlaEvidence(EvidenceTestFactory.create(), Lists.newArrayList(hla));

        LilacSummaryData lilacDataActionable = createTestLilacData("Allele 1");
        List<ProtectEvidence> evidenceActionable = hlaEvidence.evidence(lilacDataActionable);
        assertEquals(1, evidenceActionable.size());

        LilacSummaryData lilacDataNonActionable = createTestLilacData("Allele 2");
        List<ProtectEvidence> evidenceNonActionable = hlaEvidence.evidence(lilacDataNonActionable);
        assertEquals(0, evidenceNonActionable.size());
    }

    @NotNull
    private static LilacSummaryData createTestLilacData(@NotNull String hlaType) {
        List<LilacAllele> alleles = Lists.newArrayList();
        alleles.add(LilacTestFactory.builder().allele(hlaType).build());

        return ImmutableLilacSummaryData.builder().qc("PASS").alleles(alleles).build();
    }
}