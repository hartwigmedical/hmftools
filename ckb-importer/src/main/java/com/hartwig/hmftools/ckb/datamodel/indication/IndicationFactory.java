package com.hartwig.hmftools.ckb.datamodel.indication;

import com.hartwig.hmftools.ckb.json.CkbJsonDatabase;
import com.hartwig.hmftools.ckb.json.common.IndicationInfo;
import com.hartwig.hmftools.ckb.json.indication.JsonIndication;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class IndicationFactory {

    private IndicationFactory() {
    }

    @NotNull
    public static Indication extractIndication(@NotNull CkbJsonDatabase ckbJsonDatabase, @Nullable IndicationInfo indicationInfo) {
        ImmutableIndication.Builder outputBuilder = ImmutableIndication.builder();
        for (JsonIndication indication : ckbJsonDatabase.indications()) {
            if (indicationInfo.id().equals(indication.id())) {
                outputBuilder.id(indication.id())
                        .name(indication.name())
                        .source(indication.source())
                        .definition(indication.definition())
                        .currentPreferredTerm(indication.currentPreferredTerm())
                        .lastUpdateDateFromDO(indication.lastUpdateDateFromDO())
                        .altIds(indication.altIds())
                        .termId(indication.termId());
            }
        }
        return outputBuilder.build();
    }
}
