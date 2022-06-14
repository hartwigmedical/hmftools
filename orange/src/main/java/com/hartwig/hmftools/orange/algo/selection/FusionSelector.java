package com.hartwig.hmftools.orange.algo.selection;

import java.util.List;

import com.hartwig.hmftools.common.fusion.KnownFusionType;
import com.hartwig.hmftools.common.protect.ProtectEventGenerator;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.sv.linx.LinxFusion;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;

public final class FusionSelector {

    private FusionSelector() {
    }

    @NotNull
    public static List<LinxFusion> selectNonDriverFusions(@NotNull List<LinxFusion> unreportedFusions,
            @NotNull List<ProtectEvidence> evidences) {
        List<LinxFusion> filtered = Lists.newArrayList();
        for (LinxFusion fusion : unreportedFusions) {
            boolean hasEvidence = EvidenceSelector.hasEvidence(evidences, null, ProtectEventGenerator.fusionEvent(fusion));
            boolean hasReportedType = !fusion.reportedType().equals(KnownFusionType.NONE.toString());
            if (hasReportedType || hasEvidence) {
                filtered.add(fusion);
            }
        }
        return filtered;
    }
}
